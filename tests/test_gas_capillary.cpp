#include "../source/solver_lumped.h"
#include <memory>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

static void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
struct Fixture {
    using G=ConservativeGas;
    static constexpr int count=96;
    solver_lumped model{"conservation-test",""};
    solver_lumped::node inlet{},outlet{};
    solver_lumped::edge flow{},dummy{};
    std::array<std::unique_ptr<D0_edge>,8> edges;
    Fixture(bool pulmonary,double direction) {
        model.albanese_gas.enabled=true;model.conservative_gas.enabled=true;
        model.Mmax=0.;model.RQ=.7;model.edges={&dummy,&flow};
        auto& p=model.conservative_gas;
        // Closed source system: test transport and internal exchange only.
        p.k_o=0.;p.k_c=0.;
        for(auto& edge:edges) {
            edge=std::make_unique<D0_edge>("test",.001,.001,count,RBC,0.,"in","out","flow");
            edge->node_start=&inlet;edge->node_end=&outlet;edge->corr_edge=&flow;
        }
        flow.vfr=direction*.4*edges[7]->dx/.001*edges[7]->A/1.e-6;
        flow.parameter_factor=1.2;
        auto set_node=[](solver_lumped::node& n,double multiplier) {
            n.PlasmaO2_0Dn=.0025*multiplier;n.HBsat_0Dn=.8;
            n.CO2_pla_n=.024*multiplier;n.CO2_rbc_n=.018*multiplier;
            n.HCO3_pla_n=.65*multiplier;n.HCO3_rbc_n=.3*multiplier;
            n.HbCO2_n=.05*multiplier;n.RBC_fi0Dn=4.9e15*multiplier;
        };
        set_node(inlet,1.05);set_node(outlet,.95);
        model.tissueO2v.resize(count);
        for(int i=0;i<count;++i) {
            const double v=1.+.1*std::sin(i*.13);
            const G::State state{.0025*v,.8-.1*std::cos(i*.1),.024*v,.018*v,.65*v,.3*v,.05*v,.0004*v};
            for(int j=0;j<7;++j)edges[j]->fi[i]=state[j];
            edges[7]->fi[i]=4.9e15*v;model.tissueO2v[i]=state[G::T];
        }
        if(pulmonary) {
            model.pul_cap_PO2=edges[0].get();model.pul_cap_BH=edges[1].get();
            model.pul_cap_CO2_pla=edges[2].get();model.pul_cap_CO2_rbc=edges[3].get();
            model.pul_cap_HCO3_pla=edges[4].get();model.pul_cap_HCO3_rbc=edges[5].get();
            model.pul_cap_HbCO2=edges[6].get();model.pul_cap_RBC=edges[7].get();
        } else {
            model.per_cap_PO2=edges[0].get();model.per_cap_BH=edges[1].get();
            model.per_cap_CO2_pla=edges[2].get();model.per_cap_CO2_rbc=edges[3].get();
            model.per_cap_HCO3_pla=edges[4].get();model.per_cap_HCO3_rbc=edges[5].get();
            model.per_cap_HbCO2=edges[6].get();model.per_cap_RBC=edges[7].get();
        }
    }
    std::array<long double,3> cell(int i,bool tissue) const {
        const auto& p=model.conservative_gas;
        const long double o=edges[0]->fi[i]+model.Z*edges[7]->fi[i]*edges[1]->fi[i]
            +(tissue?p.phi_t/p.phi_c*model.tissueO2v[i]:0.);
        const long double c=p.phi_p()*(edges[2]->fi[i]+edges[4]->fi[i])
            +p.phi_r*(edges[3]->fi[i]+edges[5]->fi[i]+edges[6]->fi[i]);
        return {o,c,edges[7]->fi[i]};
    }
    std::array<long double,3> totals(bool tissue) const {
        std::array<long double,3> sum{};
        for(int i=0;i<count;++i) {const auto v=cell(i,tissue);for(int j=0;j<3;++j)sum[j]+=v[j];}
        return sum;
    }
};
int main() {
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    double worst_gas=0.,worst_rbc=0.;
    for(bool pulmonary:{false,true})for(double direction:{-1.,1.}) {
        Fixture serial(pulmonary,direction),parallel(pulmonary,direction);
        const auto& p=serial.model.conservative_gas;
        const auto before=serial.totals(!pulmonary);
        const auto end=serial.cell(direction>0.?Fixture::count-1:0,false);
        const auto& node=direction>0.?serial.inlet:serial.outlet;
        const std::array<long double,3> in{
            node.PlasmaO2_0Dn+serial.model.Z*node.RBC_fi0Dn*node.HBsat_0Dn,
            p.phi_p()*(node.CO2_pla_n+node.HCO3_pla_n)+p.phi_r*(node.CO2_rbc_n+node.HCO3_rbc_n+node.HbCO2_n),node.RBC_fi0Dn};
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        serial.model.advance_coupled_capillary(.001,pulmonary);
#ifdef _OPENMP
        omp_set_num_threads(6);
#endif
        parallel.model.advance_coupled_capillary(.001,pulmonary);
        for(int j=0;j<8;++j)require(serial.edges[j]->fi==parallel.edges[j]->fi,"Parallel capillary differs from serial");
        require(serial.model.tissueO2v==parallel.model.tissueO2v,"Parallel tissue differs from serial");
        const auto after=serial.totals(!pulmonary);
        for(int j=0;j<3;++j) {
            const long double error=after[j]-before[j]-.4L*(in[j]-end[j]);
            if(j<2)worst_gas=std::max(worst_gas,static_cast<double>(std::abs(error)));
            else worst_rbc=std::max(worst_rbc,static_cast<double>(std::abs(error)/before[j]));
        }
        // A worker failure must leave every capillary field uncommitted.
        parallel.edges[2]->fi[70]=-1.;
        std::array<std::vector<double>,8> snapshot;
        for(int j=0;j<8;++j)snapshot[j]=parallel.edges[j]->fi;
        const auto tissue=parallel.model.tissueO2v;
        bool rejected=false;
        const int failing_cell=direction>0.?70:69; // the injected value is advected first
        try {parallel.model.advance_coupled_capillary(.001,pulmonary);}
        catch(const std::runtime_error& ex) {
            rejected=std::string(ex.what()).find("point "+std::to_string(failing_cell))!=std::string::npos;
            if(!rejected)throw std::runtime_error(std::string("Unexpected worker diagnostic: ")+ex.what());
        }
        require(rejected,"Worker failure did not identify the failing cell");
        for(int j=0;j<8;++j)require(snapshot[j]==parallel.edges[j]->fi,"Worker failure committed a partial update");
        require(tissue==parallel.model.tissueO2v,"Worker failure committed tissue update");
    }
    require(worst_gas<1.e-8,"Capillary O2/CO2 boundary-flux balance failed");
    require(worst_rbc<1.e-14,"Capillary RBC boundary-flux balance failed");
    std::cout<<"Capillary checks: systemic/pulmonary, forward/reverse flow, varying RBC; "
             <<"serial/parallel states identical; worker rollback passed; max gas balance="
             <<worst_gas<<", RBC relative balance="<<worst_rbc<<"\n";
}
