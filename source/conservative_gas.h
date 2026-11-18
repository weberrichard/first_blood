#ifndef CONSERVATIVE_GAS_H
#define CONSERVATIVE_GAS_H

#include "albanese_gas.h"
#include <array>
#include <limits>

// Concentrations: O2 per blood volume; CO2 forms per plasma/RBC volume;
// tissue O2 per tissue volume. Hb saturation is dimensionless.
struct ConservativeGas {
    using State = std::array<double, 8>;
    enum Field { O, S, CP, CR, BP, BR, H, T };
    bool enabled = false;
    unsigned loaded = 0;
    double phi_c=.0142, phi_t=.9858, phi_r=.45;
    double beta_p=1., beta_r=.47839375648649946;
    double alpha_o_mmHg=3.11e-5, alpha_t_mmHg=3.95e-5;
    double alpha_c_water=3.27e-5, water_p=.94, water_r=.65, pa_per_mmHg=133.3224;
    double tau_o=.7469230769230769, tau_o_pul=.7469230769230769;
    double tau_c_r=.0002951690821256038, tau_b_r=.3235844356951056;
    double tau_b_reaction=.0007720610526553431, tau_h=.0703125;
    double po_alv=100., pc_alv=40., k_o=7.5006150504341366e-6, k_c=.00018183309213173663;
    double tolerance=1.e-10;
    int max_iterations=30, max_halvings=12;

    double phi_p() const { return 1.-phi_r; }
    double q_b() const { return beta_r/beta_p; }
    double q_c() const { return water_r/water_p; }
    double tau_c_p() const { return tau_c_r*phi_p()/(phi_r*q_c()); }
    double tau_b_p() const { return tau_b_r*phi_p()/(phi_r*q_b()); }
    double alpha_c_p(double vm) const { return alpha_c_water*vm*water_p/pa_per_mmHg; }
    double alpha_c_r(double vm) const { return alpha_c_water*vm*water_r/pa_per_mmHg; }

    bool read_row(const std::vector<std::string>& row) {
        if(row.empty()) return false;
        const std::string& key=row[0];
        auto values=[&](std::initializer_list<double*> targets, unsigned bit) {
            if(row.size()!=targets.size()+1) throw std::runtime_error(key+": incorrect parameter count");
            std::size_t i=1;
            for(double* target:targets) *target=AlbaneseGas::number(row[i++]);
            loaded|=bit;
        };
        if(key=="gas_solver") {
            if(row.size()!=2 || (row[1]!="implicit" && row[1]!="legacy"))
                throw std::runtime_error("gas_solver requires implicit|legacy");
            enabled=row[1]=="implicit";
        } else if(key=="gas_fractions") values({&phi_c,&phi_t,&phi_r},1);
        else if(key=="gas_partition") values({&beta_p,&beta_r},2);
        else if(key=="gas_solubility") values({&alpha_o_mmHg,&alpha_t_mmHg,&alpha_c_water,&water_p,&water_r,&pa_per_mmHg},4);
        else if(key=="gas_relaxation") values({&tau_o,&tau_o_pul,&tau_c_r,&tau_b_r,&tau_b_reaction,&tau_h},8);
        else if(key=="gas_pulmonary") values({&po_alv,&pc_alv,&k_o,&k_c},16);
        else if(key=="gas_numerics") {
            double iterations=0.,halvings=0.;
            values({&tolerance,&iterations,&halvings},32);
            if(iterations!=std::floor(iterations) || iterations<1 || iterations>100 ||
               halvings!=std::floor(halvings) || halvings<0 || halvings>20)
                throw std::runtime_error("gas_numerics: invalid iteration/subdivision limit");
            max_iterations=static_cast<int>(iterations); max_halvings=static_cast<int>(halvings);
        } else return false;
        return true;
    }

    void validate() const {
        if(!enabled) return;
        if((loaded&31)!=31) throw std::runtime_error("Implicit gas solver requires gas_fractions, gas_partition, gas_solubility, gas_relaxation and gas_pulmonary rows");
        for(double x:{phi_c,phi_t,phi_r,beta_p,beta_r,alpha_o_mmHg,alpha_t_mmHg,
            alpha_c_water,water_p,water_r,pa_per_mmHg,tau_o,tau_o_pul,tau_c_r,
            tau_b_r,tau_b_reaction,tau_h,po_alv,pc_alv,k_o,k_c,tolerance})
            if(!std::isfinite(x)||x<=0.) throw std::runtime_error("Implicit gas parameters must be finite and positive");
        if(phi_r>=1. || phi_c>=1. || water_p>1. || water_r>1. ||
           std::abs(phi_c+phi_t-1.)>1.e-10 || tolerance>1.e-4)
            throw std::runtime_error("Implicit gas fractions or tolerance are invalid");
    }

    struct Environment {
        bool pulmonary=false;
        double capacity=0., mmax=0., c50=0., rq=0.;
        double oxygen_conductance=0., profile=0.; // conductance in 1/(Pa s)
    };

    State rates(const State& y, const Environment& e, const AlbaneseGas& gas) const {
        const double po=y[O]/alpha_o_mmHg;
        const double ap=alpha_c_p(gas.molar_volume_l_mol), ar=alpha_c_r(gas.molar_volume_l_mol);
        const double pc=y[CP]/ap;
        const double equilibrium=gas.saturation(y[O],po,pc/pa_per_mmHg,e.capacity);
        const double departure=y[S]-equilibrium;
        const double tau=e.pulmonary?tau_o_pul:tau_o;
        const double jo=e.capacity*departure/tau;
        const auto target=gas.carbon_targets(po,pc,ap,ar,phi_p(),phi_r,pa_per_mmHg);
        const double jc=phi_r/tau_c_r*(q_c()*y[CP]-y[CR]);
        const double jb=phi_r/tau_b_r*(q_b()*y[BP]-y[BR]);
        const double rb=phi_r/tau_b_reaction*(target.bicarbonate-phi_p()*y[BP]-phi_r*y[BR]);
        const double rh=phi_r/tau_h*(target.carbamino-phi_r*y[H]);
        State f{};
        f[O]=jo; f[S]=-departure/tau;
        f[CP]=-jc/phi_p(); f[CR]=(jc-rb-rh)/phi_r;
        f[BP]=-jb/phi_p(); f[BR]=(jb+rb)/phi_r; f[H]=rh/phi_r;
        if(e.pulmonary) {
            f[O]+=k_o*e.profile*(po_alv-po)*pa_per_mmHg;
            f[CP]+=k_c*e.profile*(pc_alv*pa_per_mmHg-pc);
        } else {
            const double diffusion=e.oxygen_conductance*pa_per_mmHg*(po-y[T]/alpha_t_mmHg);
            const double metabolism=e.mmax*y[T]/(y[T]+e.c50);
            f[O]-=diffusion;
            f[T]=phi_c/phi_t*diffusion-metabolism;
            // Exactly the same metabolism at the same iterate feeds both gases.
            f[CP]+=e.rq*phi_t/(phi_c*phi_p())*metabolism;
        }
        return f;
    }

    static bool physical(const State& y) {
        for(double x:y) if(!std::isfinite(x)||x<0.) return false;
        return y[S]<=1.;
    }

    bool backward_euler(State& y, double dt, const Environment& e, const AlbaneseGas& gas) const {
        const State old=y;
        auto residual=[&](const State& x) {
            const State f=rates(x,e,gas); State r{};
            for(int j=0;j<8;++j) r[j]=x[j]-old[j]-dt*f[j];
            return r;
        };
        auto norm=[&](const State& r,const State& x) {
            double n=0.;
            for(int j=0;j<8;++j) n=std::max(n,std::abs(r[j])/std::max({1.e-4,std::abs(old[j]),std::abs(x[j])}));
            return n;
        };
        try {
            for(int iteration=0;iteration<max_iterations;++iteration) {
                const State r=residual(y); const double n=norm(r,y);
                if(n<tolerance) return true;
                std::array<std::array<double,9>,8> a{};
                for(int j=0;j<8;++j) {
                    State x=y; double h=1.e-7*std::max(1.e-3,std::abs(y[j]));
                    if(j==S && x[j]+h>1.) h=-h;
                    x[j]+=h; const State shifted=residual(x);
                    for(int i=0;i<8;++i) a[i][j]=(shifted[i]-r[i])/h;
                }
                for(int i=0;i<8;++i) a[i][8]=-r[i];
                for(int j=0;j<8;++j) {
                    int pivot=j;
                    for(int i=j+1;i<8;++i) if(std::abs(a[i][j])>std::abs(a[pivot][j])) pivot=i;
                    if(std::abs(a[pivot][j])<1.e-18) return false;
                    std::swap(a[j],a[pivot]);
                    const double d=a[j][j]; for(int k=j;k<9;++k) a[j][k]/=d;
                    for(int i=0;i<8;++i) if(i!=j) {
                        const double factor=a[i][j];
                        for(int k=j;k<9;++k) a[i][k]-=factor*a[j][k];
                    }
                }
                bool accepted=false;
                for(double damping=1.;damping>=1./65536.;damping*=.5) {
                    State trial=y;
                    for(int j=0;j<8;++j) trial[j]+=damping*a[j][8];
                    if(!physical(trial)) continue;
                    try {
                        if(norm(residual(trial),trial)<n) { y=trial; accepted=true; break; }
                    } catch(const std::runtime_error&) {}
                }
                if(!accepted) return false;
            }
        } catch(const std::runtime_error&) { return false; }
        return false;
    }

    void advance(State& y,double dt,const Environment& e,const AlbaneseGas& gas,int depth=0) const {
        if(!std::isfinite(dt)||dt<0. || !physical(y)) throw std::runtime_error("Invalid implicit gas state/timestep");
        if(dt==0.) return;
        State trial=y;
        if(backward_euler(trial,dt,e,gas)) {y=trial;return;}
        if(depth>=max_halvings) throw std::runtime_error("Implicit coupled gas solve failed to converge");
        // Retry a failed step without committing a partial update.
        trial=y;
        advance(trial,dt*.5,e,gas,depth+1);
        advance(trial,dt*.5,e,gas,depth+1);
        y=trial;
    }
};
#endif
