#include "../source/solver_lumped.h"
#include <memory>
#include <iostream>
#include <stdexcept>
static void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(){
    using N=solver_lumped::node;
    struct Case {TransportType type;std::vector<D0_edge*> N::* incoming;std::vector<D0_edge*> N::* outgoing;double N::* value;};
    const Case cases[]={
        {RBC,&N::D0_edges_in_RBC,&N::D0_edges_out_RBC,&N::RBC_fi0Dn},
        {HB_O2_saturation,&N::D0_edges_in_HBsat,&N::D0_edges_out_HBsat,&N::HBsat_0Dn},
        {C_Plasma_O2,&N::D0_edges_in_PlasmaO2,&N::D0_edges_out_PlasmaO2,&N::PlasmaO2_0Dn},
        {CO2_pla,&N::D0_edges_in_CO2_pla,&N::D0_edges_out_CO2_pla,&N::CO2_pla_n},
        {CO2_rbc,&N::D0_edges_in_CO2_rbc,&N::D0_edges_out_CO2_rbc,&N::CO2_rbc_n},
        {HCO3_pla,&N::D0_edges_in_HCO3_pla,&N::D0_edges_out_HCO3_pla,&N::HCO3_pla_n},
        {HCO3_rbc,&N::D0_edges_in_HCO3_rbc,&N::D0_edges_out_HCO3_rbc,&N::HCO3_rbc_n},
        {HbCO2,&N::D0_edges_in_HbCO2,&N::D0_edges_out_HbCO2,&N::HbCO2_n}};
    for(const auto& test:cases) {
        solver_lumped model("mixing-test","");N node{};model.nodes={&node};
        for(const auto& item:cases)node.*item.value=-9.;
        solver_lumped::edge positive{},negative{},closed{};
        positive.vfr=2.;negative.vfr=-1.;closed.vfr=100.;closed.is_open=false;
        D0_edge in("in",1.,1.,3,test.type,.8,"","",""),out("out",1.,1.,3,test.type,.2,"","",""),blocked("closed",1.,1.,3,test.type,99.,"","","");
        in.corr_edge=&positive;out.corr_edge=&negative;blocked.corr_edge=&closed;
        node.*test.incoming={&in,&blocked};node.*test.outgoing={&out};
        D0_transport transport(test.type);transport.update_nodes(model);
        require(std::abs(node.*test.value-.6)<1.e-15,"Typed mixing changed the flow-weighted concentration");
        require(std::abs(3.*(node.*test.value)-(2.*.8+.2))<1.e-14,"Mixed species inflow balance failed");
        for(const auto& other:cases)if(other.type!=test.type)require(node.*other.value==-9.,"Typed mixing wrote a different species");
        node.*test.value=.33;node.is_master_node=true;transport.update_nodes(model);
        require(node.*test.value==.33,"Master node was overwritten");
        node.is_master_node=false;positive.vfr=-2.;negative.vfr=1.;transport.update_nodes(model);
        require(node.*test.value==.33,"Zero incoming flow did not retain concentration");
        positive.vfr=0.;negative.vfr=0.;transport.update_nodes(model);
        require(node.*test.value==.33,"Stagnant node changed concentration");
    }
    std::cout<<"Typed mixing: all eight species, inflow conservation, closed edges, reversed/stagnant flow and master-node checks passed\n";
}
