#include "../source/conservative_gas.h"
#include <iostream>
#include <random>
#ifdef GAS_REFERENCE_HEADER
#undef CONSERVATIVE_GAS_H
#define ConservativeGas ReferenceConservativeGas
#include GAS_REFERENCE_HEADER
#undef ConservativeGas
#endif

static void require(bool ok,const char* why) { if(!ok)throw std::runtime_error(why); }
static double carbon(const ConservativeGas::State& y,const ConservativeGas& p) {
    return p.phi_p()*(y[p.CP]+y[p.BP])+p.phi_r*(y[p.CR]+y[p.BR]+y[p.H]);
}
static double oxygen(const ConservativeGas::State& y,const ConservativeGas& p,double capacity) {
    return y[p.O]+capacity*y[p.S]+p.phi_t/p.phi_c*y[p.T];
}
int main() {
    using G=ConservativeGas;
    G p;AlbaneseGas gas;gas.enabled=true;
    std::mt19937 random(20261119);std::uniform_real_distribution<double> uniform(0.,1.);
    double worst_derivative=0.,worst_balance=0.,worst_reference=0.,worst_residual=0.;
    int comparisons=0,steps=0,reduced_steps=0;
    double worst_reduced_derivative=0.;
    for(int sample=0;sample<160;++sample) {
        G::Environment e;e.pulmonary=sample%2;e.capacity=.15+.1*uniform(random);
        e.mmax=1.e-4+2.e-4*uniform(random);e.c50=2.6e-5;e.rq=.7;
        e.oxygen_conductance=1.e-4*uniform(random);e.profile=uniform(random);
        const double po=10.+140.*uniform(random),pc=10.+60.*uniform(random);
        G::State y{po*p.alpha_o_mmHg,.05+.9*uniform(random),pc*p.pa_per_mmHg*p.alpha_c_p(gas.molar_volume_l_mol),0.,0.,0.,0.,(.5+60.*uniform(random))*p.alpha_t_mmHg};
        const auto target=gas.carbon_targets(po,pc*p.pa_per_mmHg,p.alpha_c_p(gas.molar_volume_l_mol),p.alpha_c_r(gas.molar_volume_l_mol),p.phi_p(),p.phi_r,p.pa_per_mmHg);
        y[G::CR]=y[G::CP]*p.q_c()*(.8+.4*uniform(random));
        y[G::BP]=target.bicarbonate/(p.phi_p()+p.phi_r*p.q_b())*(.8+.4*uniform(random));
        y[G::BR]=p.q_b()*y[G::BP]*(.8+.4*uniform(random));
        y[G::H]=target.carbamino/p.phi_r*(.8+.4*uniform(random));
        G::Jacobian jac;p.rates(y,e,gas,&jac);
        for(int j=0;j<8;++j) {
            const double h=std::max(1.e-9,std::abs(y[j])*1.e-5);
            auto plus=y,minus=y;plus[j]+=h;minus[j]-=h;
            const auto fp=p.rates(plus,e,gas),fm=p.rates(minus,e,gas);
            for(int i=0;i<8;++i) {
                const double numeric=(fp[i]-fm[i])/(2.*h);
                worst_derivative=std::max(worst_derivative,std::abs(numeric-jac[i][j])/std::max({1.,std::abs(numeric),std::abs(jac[i][j])}));
            }
        }
        for(double dt:{1.e-6,1.e-4,1.e-3,5.e-3}) {
            auto next=y;
            const G::Reduced z{y[G::O],y[G::CP],y[G::T]};
            const auto reduced=p.reduced_evaluation(y,z,dt,e,gas,true);
            for(int j=0;j<3;++j) {
                const double h=std::max(1.e-9,std::abs(z[j])*1.e-5);
                auto plus=z,minus=z;plus[j]+=h;minus[j]-=h;
                const auto fp=p.reduced_evaluation(y,plus,dt,e,gas,false);
                const auto fm=p.reduced_evaluation(y,minus,dt,e,gas,false);
                for(int i=0;i<3;++i) {
                    const double numeric=(fp.residual[i]-fm.residual[i])/(2.*h);
                    worst_reduced_derivative=std::max(worst_reduced_derivative,std::abs(numeric-reduced.jacobian[i][j])/std::max({1.,std::abs(numeric),std::abs(reduced.jacobian[i][j])}));
                }
            }
            auto direct_reduced=y;
            if(p.backward_euler_reduced(direct_reduced,dt,e,gas))++reduced_steps;
            // The comparison is against the original full backward-Euler
            // equations, including a direct check of every residual component.
            const bool direct=p.backward_euler(next,dt,e,gas);
            if(!direct) {next=y;p.advance(next,dt,e,gas);}
            require(G::physical(next),"Optimization lost physical state bounds");
            if(direct) {
                const auto rate=p.rates(next,e,gas);
                for(int j=0;j<8;++j)worst_residual=std::max(worst_residual,std::abs(next[j]-y[j]-dt*rate[j])/std::max({1.e-4,std::abs(next[j]),std::abs(y[j])}));
                const double metabolism=e.mmax*next[G::T]/(next[G::T]+e.c50);
                const double flux_o=e.pulmonary?p.k_o*e.profile*(p.po_alv-next[G::O]/p.alpha_o_mmHg)*p.pa_per_mmHg:-p.phi_t/p.phi_c*metabolism;
                const double flux_c=e.pulmonary?p.phi_p()*p.k_c*e.profile*(p.pc_alv*p.pa_per_mmHg-next[G::CP]/p.alpha_c_p(gas.molar_volume_l_mol)):e.rq*p.phi_t/p.phi_c*metabolism;
                worst_balance=std::max({worst_balance,std::abs(oxygen(next,p,e.capacity)-oxygen(y,p,e.capacity)-dt*flux_o),std::abs(carbon(next,p)-carbon(y,p)-dt*flux_c)});
            }
#ifdef GAS_REFERENCE_HEADER
            ReferenceConservativeGas reference;
            ReferenceConservativeGas::Environment re;
            re.pulmonary=e.pulmonary;re.capacity=e.capacity;re.mmax=e.mmax;re.c50=e.c50;
            re.rq=e.rq;re.oxygen_conductance=e.oxygen_conductance;re.profile=e.profile;
            auto expected=y;reference.advance(expected,dt,re,gas);
            for(int j=0;j<8;++j)worst_reference=std::max(worst_reference,std::abs(next[j]-expected[j])/std::max({1.e-4,std::abs(next[j]),std::abs(expected[j])}));
            ++comparisons;
#endif
            ++steps;
        }
    }
    require(worst_derivative<5.e-6,"Analytic Jacobian differs from central differences");
    require(worst_balance<1.e-9,"Optimization violates integrated O2/CO2 balance");
    require(worst_residual<2.e-10,"Optimization violates original backward-Euler equations");
    require(worst_reference<1.e-6,"Optimization differs from reference source solution");
    require(reduced_steps>steps*.9,"Reduced solve rarely converges without full fallback");
    require(worst_reduced_derivative<5.e-6,"Reduced Jacobian differs from central differences");
    // Zero-pressure stationary states use the full fallback rather than an
    // undefined derivative of the empirical carbon curve at zero pressure.
    G::Environment closed;closed.capacity=.2;closed.c50=2.6e-5;
    G::State empty{};p.advance(empty,.001,closed,gas);
    for(double value:empty)require(value==0.,"Empty closed system changed");
    G::State no_carbon{};no_carbon[G::O]=.002;
    no_carbon[G::S]=gas.saturation(no_carbon[G::O],no_carbon[G::O]/p.alpha_o_mmHg,0.,closed.capacity);
    const auto stationary=no_carbon;p.advance(no_carbon,.001,closed,gas);
    for(int i=0;i<8;++i)require(std::abs(no_carbon[i]-stationary[i])<1.e-12,"Zero-carbon equilibrium changed");
    std::cout<<"Optimization checks: "<<steps<<" source steps, "<<comparisons
             <<" reference comparisons; max Jacobian error="<<worst_derivative
             <<", O2/CO2 balance error="<<worst_balance
             <<", original equation residual="<<worst_residual
             <<", reference difference="<<worst_reference
             <<", direct reduced solves="<<reduced_steps
             <<", reduced Jacobian error="<<worst_reduced_derivative<<"\n";
}
