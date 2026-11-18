#include "../source/conservative_gas.h"
#include <iostream>

static void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
static ConservativeGas parameters() {
    ConservativeGas p;
    for(const auto& row:std::vector<std::vector<std::string>>{
        {"gas_solver","implicit"},{"gas_fractions","0.0142","0.9858","0.45"},
        {"gas_partition","1","0.47839375648649946"},
        {"gas_solubility","3.11e-5","3.95e-5","3.27e-5",".94",".65","133.3224"},
        {"gas_relaxation",".7469230769230769",".7469230769230769",".0002951690821256038",".3235844356951056",".0007720610526553431",".0703125"},
        {"gas_pulmonary","100","40","7.5006150504341366e-6",".00018183309213173663"}})
        require(p.read_row(row),"New parameter row not recognized");
    p.validate();return p;
}
static double carbon(const ConservativeGas::State& y,const ConservativeGas& p) {
    return p.phi_p()*(y[p.CP]+y[p.BP])+p.phi_r*(y[p.CR]+y[p.BR]+y[p.H]);
}
static double oxygen(const ConservativeGas::State& y,const ConservativeGas& p,double capacity) {
    return y[p.O]+capacity*y[p.S]+p.phi_t/p.phi_c*y[p.T];
}
int main() {
    using G=ConservativeGas;
    G p=parameters();AlbaneseGas gas;gas.enabled=true;
    G::Environment e;e.capacity=4.9e15*3.73e-17;e.c50=2.6e-5;e.rq=.7;
    G::State initial{.003,.92,.026,.04,.7,.09,.04,.0004};
    const double c0=carbon(initial,p),o0=oxygen(initial,p,e.capacity);
    auto y=initial;
    // Ten milliseconds, with a 1 ms step larger than the fast exchange time.
    for(int i=0;i<10;++i)p.advance(y,.001,e,gas);
    require(G::physical(y),"Stiff source solve lost positivity/saturation bounds");
    require(std::abs(carbon(y,p)-c0)<1.e-9,"Closed source solve created carbon");
    require(std::abs(oxygen(y,p,e.capacity)-o0)<1.e-9,"Closed source solve created oxygen");
    require(std::abs(y[p.CP]-initial[p.CP])>1.e-5,"Source solve did not evolve");
    // A local equilibrium must remain stationary, including all carbon forms.
    G::State eq=initial;
    eq[G::O]=60.*p.alpha_o_mmHg;
    eq[G::S]=gas.saturation(eq[G::O],60.,40.,e.capacity);
    eq[G::CP]=40.*p.pa_per_mmHg*p.alpha_c_p(gas.molar_volume_l_mol);
    eq[G::CR]=p.q_c()*eq[G::CP];
    const auto target=gas.carbon_targets(60.,40.*p.pa_per_mmHg,p.alpha_c_p(gas.molar_volume_l_mol),p.alpha_c_r(gas.molar_volume_l_mol),p.phi_p(),p.phi_r,p.pa_per_mmHg);
    eq[G::BP]=target.bicarbonate/(p.phi_p()+p.phi_r*p.q_b());
    eq[G::BR]=p.q_b()*eq[G::BP];eq[G::H]=target.carbamino/p.phi_r;
    auto stationary=eq;p.advance(stationary,.005,e,gas);
    for(int i=0;i<8;++i)require(std::abs(stationary[i]-eq[i])<1.e-10,"Equilibrium drifted");
    // Independent unit-invariance test: common beta scaling must change nothing.
    G scaled=p;scaled.beta_p*=.01;scaled.beta_r*=.01;
    auto a=initial,b=initial;p.advance(a,.001,e,gas);scaled.advance(b,.001,e,gas);
    for(int i=0;i<8;++i)require(std::abs(a[i]-b[i])<1.e-10,"Common beta scale changed physics");
    require(std::abs(p.phi_r*p.q_c()/p.tau_c_r-p.phi_p()/p.tau_c_p())<1.e-9,"CO2 partner time violates mass balance");
    require(std::abs(p.phi_r*p.q_b()/p.tau_b_r-p.phi_p()/p.tau_b_p())<1.e-12,"Bicarbonate partner time violates mass balance");
    // Metabolism must consume O2 and create precisely RQ times as much CO2.
    e.mmax=1.4e-4;e.oxygen_conductance=1.e-4;
    y=initial;p.advance(y,.001,e,gas);
    const double consumed=o0-oxygen(y,p,e.capacity);
    require(consumed>0.,"Metabolism did not consume oxygen");
    require(std::abs((carbon(y,p)-c0)-e.rq*consumed)<1.e-9,"O2/CO2 metabolic stoichiometry is inconsistent");
    // Backward Euler pulmonary balance against independently integrated external flux.
    e.pulmonary=true;e.profile=.7;y=initial;p.advance(y,.001,e,gas);
    const double po=y[p.O]/p.alpha_o_mmHg,pc=y[p.CP]/p.alpha_c_p(gas.molar_volume_l_mol);
    require(std::abs(oxygen(y,p,e.capacity)-o0-.001*p.k_o*e.profile*(p.po_alv-po)*p.pa_per_mmHg)<1.e-9,"Pulmonary oxygen balance mismatch");
    require(std::abs(carbon(y,p)-c0-.001*p.phi_p()*p.k_c*e.profile*(p.pc_alv*p.pa_per_mmHg-pc))<1.e-9,"Pulmonary carbon balance mismatch");
    for(auto row:std::vector<std::vector<std::string>>{
        {"gas_solver","unknown"},{"gas_partition","1"},{"gas_numerics","1e-10","1.5","12"},
        {"gas_relaxation","nan","1","1","1","1","1"}}) {
        bool rejected=false;try{G bad;bad.read_row(row);}catch(const std::exception&){rejected=true;}
        require(rejected,"Malformed file row was accepted");
    }
    bool rejected=false;try{G bad;bad.enabled=true;bad.validate();}catch(const std::exception&){rejected=true;}
    require(rejected,"Incomplete implicit parameter set was accepted");
    rejected=false;try{G bad=p;bad.phi_r=1.;bad.validate();}catch(const std::exception&){rejected=true;}
    require(rejected,"Invalid fraction was accepted");
    std::cout<<"Implicit gas: parsing, stiff relaxation, positivity, beta scaling, conservation and metabolism passed\n";
}
