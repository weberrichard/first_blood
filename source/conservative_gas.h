#ifndef CONSERVATIVE_GAS_H
#define CONSERVATIVE_GAS_H

#include "albanese_gas.h"
#include <array>
#include <limits>

// Concentrations: O2 per blood volume; CO2 forms per plasma/RBC volume;
// tissue O2 per tissue volume. Hb saturation is dimensionless.
struct ConservativeGas {
    using State = std::array<double, 8>;
    using Jacobian = std::array<State, 8>;
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

    State rates(const State& y, const Environment& e, const AlbaneseGas& gas,
                Jacobian* jacobian=nullptr) const {
        const double po=y[O]/alpha_o_mmHg;
        const double ap=alpha_c_p(gas.molar_volume_l_mol), ar=alpha_c_r(gas.molar_volume_l_mol);
        const double pc=y[CP]/ap;
        AlbaneseGas::OxygenGradient og;
        AlbaneseGas::CarbonGradient cg;
        const double equilibrium=gas.saturation(y[O],po,pc/pa_per_mmHg,e.capacity,jacobian?&og:nullptr);
        const double departure=y[S]-equilibrium;
        const double tau=e.pulmonary?tau_o_pul:tau_o;
        const double jo=e.capacity*departure/tau;
        const auto target=gas.carbon_targets(po,pc,ap,ar,phi_p(),phi_r,pa_per_mmHg,jacobian?&cg:nullptr);
        const double jc=phi_r/tau_c_r*(q_c()*y[CP]-y[CR]);
        const double jb=phi_r/tau_b_r*(q_b()*y[BP]-y[BR]);
        const double rb=phi_r/tau_b_reaction*(target.bicarbonate-phi_p()*y[BP]-phi_r*y[BR]);
        const double rh=phi_r/tau_h*(target.carbamino-phi_r*y[H]);
        State f{};
        f[O]=jo; f[S]=-departure/tau;
        f[CP]=-jc/phi_p(); f[CR]=(jc-rb-rh)/phi_r;
        f[BP]=-jb/phi_p(); f[BR]=(jb+rb)/phi_r; f[H]=rh/phi_r;
        if(jacobian) {
            *jacobian=Jacobian{};
            State ds{},dc{},db{},drb{},drh{};
            ds[O]=-og.pressure/alpha_o_mmHg-og.plasma;
            ds[CP]=-og.other/(ap*pa_per_mmHg);ds[S]=1.;
            dc[CP]=phi_r/tau_c_r*q_c();dc[CR]=-phi_r/tau_c_r;
            db[BP]=phi_r/tau_b_r*q_b();db[BR]=-phi_r/tau_b_r;
            drb[O]=phi_r/tau_b_reaction*cg.bicarbonate_po/alpha_o_mmHg;
            drb[CP]=phi_r/tau_b_reaction*cg.bicarbonate_pc/ap;
            drb[BP]=-phi_r/tau_b_reaction*phi_p();drb[BR]=-phi_r/tau_b_reaction*phi_r;
            drh[O]=phi_r/tau_h*cg.carbamino_po/alpha_o_mmHg;
            drh[CP]=phi_r/tau_h*cg.carbamino_pc/ap;drh[H]=-phi_r*phi_r/tau_h;
            for(int j=0;j<8;++j) {
                (*jacobian)[O][j]=e.capacity/tau*ds[j];(*jacobian)[S][j]=-ds[j]/tau;
                (*jacobian)[CP][j]=-dc[j]/phi_p();(*jacobian)[CR][j]=(dc[j]-drb[j]-drh[j])/phi_r;
                (*jacobian)[BP][j]=-db[j]/phi_p();(*jacobian)[BR][j]=(db[j]+drb[j])/phi_r;
                (*jacobian)[H][j]=drh[j]/phi_r;
            }
        }
        if(e.pulmonary) {
            f[O]+=k_o*e.profile*(po_alv-po)*pa_per_mmHg;
            f[CP]+=k_c*e.profile*(pc_alv*pa_per_mmHg-pc);
            if(jacobian) {
                (*jacobian)[O][O]-=k_o*e.profile*pa_per_mmHg/alpha_o_mmHg;
                (*jacobian)[CP][CP]-=k_c*e.profile/ap;
            }
        } else {
            const double diffusion=e.oxygen_conductance*pa_per_mmHg*(po-y[T]/alpha_t_mmHg);
            const double metabolism=e.mmax*y[T]/(y[T]+e.c50);
            f[O]-=diffusion;
            f[T]=phi_c/phi_t*diffusion-metabolism;
            // Exactly the same metabolism at the same iterate feeds both gases.
            f[CP]+=e.rq*phi_t/(phi_c*phi_p())*metabolism;
            if(jacobian) {
                const double d_o=e.oxygen_conductance*pa_per_mmHg/alpha_o_mmHg;
                const double d_t=-e.oxygen_conductance*pa_per_mmHg/alpha_t_mmHg;
                const double d_m=e.mmax*e.c50/((y[T]+e.c50)*(y[T]+e.c50));
                (*jacobian)[O][O]-=d_o;(*jacobian)[O][T]-=d_t;
                (*jacobian)[T][O]=phi_c/phi_t*d_o;(*jacobian)[T][T]=phi_c/phi_t*d_t-d_m;
                (*jacobian)[CP][T]+=e.rq*phi_t/(phi_c*phi_p())*d_m;
            }
        }
        return f;
    }

    static bool physical(const State& y) {
        for(double x:y) if(!std::isfinite(x)||x<0.) return false;
        return y[S]<=1.;
    }

    bool backward_euler_full(State& y, double dt, const Environment& e, const AlbaneseGas& gas) const {
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
                const bool analytic=y[O]>0. && y[CP]>0.;
                Jacobian jacobian{};
                const State f=rates(y,e,gas,analytic?&jacobian:nullptr);
                State r{};for(int i=0;i<8;++i)r[i]=y[i]-old[i]-dt*f[i];
                const double n=norm(r,y);
                if(n<tolerance) return true;
                std::array<std::array<double,9>,8> a{};
                if(analytic) {
                    for(int i=0;i<8;++i)for(int j=0;j<8;++j)
                        a[i][j]=(i==j?1.:0.)-dt*jacobian[i][j];
                } else for(int j=0;j<8;++j) {
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

    using Reduced = std::array<double,3>; // O, Cp, stationary tissue O2
    struct ReducedEvaluation {
        State state{};
        Reduced residual{};
        std::array<Reduced,3> jacobian{};
    };

    ReducedEvaluation reduced_evaluation(const State& old,const Reduced& z,double dt,
                                         const Environment& e,const AlbaneseGas& gas,
                                         bool derivatives) const {
        ReducedEvaluation out;
        auto& y=out.state;y=old;y[O]=z[0];y[CP]=z[1];y[T]=z[2];
        const double ap=alpha_c_p(gas.molar_volume_l_mol),ar=alpha_c_r(gas.molar_volume_l_mol);
        const double po=y[O]/alpha_o_mmHg,pc=y[CP]/ap;
        AlbaneseGas::OxygenGradient og;
        AlbaneseGas::CarbonGradient cg;
        const double seq=gas.saturation(y[O],po,pc/pa_per_mmHg,e.capacity,derivatives?&og:nullptr);
        const auto target=gas.carbon_targets(po,pc,ap,ar,phi_p(),phi_r,pa_per_mmHg,derivatives?&cg:nullptr);
        const double tau=e.pulmonary?tau_o_pul:tau_o;
        const double sgain=dt/(tau+dt),hgain=dt/(tau_h+dt*phi_r);
        y[S]=(tau*old[S]+dt*seq)/(tau+dt);
        y[H]=(tau_h*old[H]+dt*target.carbamino)/(tau_h+dt*phi_r);

        // Exact elimination of the two linear bicarbonate BE equations.
        // The expanded determinant has positive terms, avoiding cancellation.
        const double a=dt*phi_r/(phi_p()*tau_b_r),b=dt/tau_b_r,c=dt/tau_b_reaction;
        const double a11=1.+a*q_b(),a21=-b*q_b()+c*phi_p(),a22=1.+b+c*phi_r;
        const double determinant=1.+b+a*q_b()+c*(phi_r+a*(q_b()*phi_r+phi_p()));
        const double rhs=old[BR]+c*target.bicarbonate;
        y[BP]=(a22*old[BP]+a*rhs)/determinant;
        y[BR]=(a11*rhs-a21*old[BP])/determinant;

        // The RBC free-CO2 equation is also linear once O and Cp are given.
        const double d=dt/tau_c_r;
        y[CR]=(old[CR]+d*q_c()*y[CP]
            -c*(target.bicarbonate-phi_p()*y[BP]-phi_r*y[BR])
            -hgain*(target.carbamino-phi_r*old[H]))/(1.+d);
        const double exchange=dt*phi_r/(phi_p()*tau_c_r);
        out.residual={y[O]-old[O]+e.capacity*(y[S]-old[S]),
                      y[CP]-old[CP]+exchange*(q_c()*y[CP]-y[CR]),y[T]-old[T]};
        if(derivatives) {
            const double so=sgain*(og.pressure/alpha_o_mmHg+og.plasma);
            const double sc=sgain*og.other/(ap*pa_per_mmHg);
            const double bgain=c*(1.+b+a*q_b())/determinant;
            const double ro=(-bgain*cg.bicarbonate_po-hgain*cg.carbamino_po)/(alpha_o_mmHg*(1.+d));
            const double rc=(d*q_c()-(bgain*cg.bicarbonate_pc+hgain*cg.carbamino_pc)/ap)/(1.+d);
            out.jacobian={Reduced{1.+e.capacity*so,e.capacity*sc,0.},
                          Reduced{-exchange*ro,1.+exchange*(q_c()-rc),0.},Reduced{0.,0.,1.}};
        }
        if(e.pulmonary) {
            out.residual[0]-=dt*k_o*e.profile*(po_alv-po)*pa_per_mmHg;
            out.residual[1]-=dt*k_c*e.profile*(pc_alv*pa_per_mmHg-pc);
            if(derivatives) {
                out.jacobian[0][0]+=dt*k_o*e.profile*pa_per_mmHg/alpha_o_mmHg;
                out.jacobian[1][1]+=dt*k_c*e.profile/ap;
            }
        } else {
            const double conductance=e.oxygen_conductance*pa_per_mmHg;
            const double diffusion=conductance*(po-y[T]/alpha_t_mmHg);
            const double metabolism=e.mmax*y[T]/(y[T]+e.c50);
            out.residual[0]+=dt*diffusion;
            out.residual[1]-=dt*e.rq*phi_t/(phi_c*phi_p())*metabolism;
            out.residual[2]+=dt*(metabolism-phi_c/phi_t*diffusion);
            if(derivatives) {
                const double dm=e.mmax*e.c50/((y[T]+e.c50)*(y[T]+e.c50));
                out.jacobian[0][0]+=dt*conductance/alpha_o_mmHg;
                out.jacobian[0][2]-=dt*conductance/alpha_t_mmHg;
                out.jacobian[1][2]-=dt*e.rq*phi_t/(phi_c*phi_p())*dm;
                out.jacobian[2][0]-=dt*phi_c/phi_t*conductance/alpha_o_mmHg;
                out.jacobian[2][2]+=dt*(dm+phi_c/phi_t*conductance/alpha_t_mmHg);
            }
        }
        return out;
    }

    bool backward_euler_reduced(State& y,double dt,const Environment& e,const AlbaneseGas& gas) const {
        const State old=y;
        if(y[O]<=0. || y[CP]<=0.)return false; // singular curve derivative at zero pressure
        Reduced z{y[O],y[CP],y[T]};
        const int fields[3]={O,CP,T};
        auto norm=[&](const ReducedEvaluation& value,const Reduced& x) {
            double n=0.;
            for(int i=0;i<3;++i)n=std::max(n,std::abs(value.residual[i])/
                std::max({1.e-4,std::abs(old[fields[i]]),std::abs(x[i])}));
            return n;
        };
        try {
            for(int iteration=0;iteration<max_iterations;++iteration) {
                const auto value=reduced_evaluation(old,z,dt,e,gas,true);
                const double n=norm(value,z);
                if(n<tolerance) {
                    if(!physical(value.state))return false;
                    // Verify ALL eight original equations before committing.
                    const auto f=rates(value.state,e,gas);
                    for(int i=0;i<8;++i)
                        if(std::abs(value.state[i]-old[i]-dt*f[i])/
                           std::max({1.e-4,std::abs(old[i]),std::abs(value.state[i])})>=tolerance)return false;
                    y=value.state;return true;
                }
                std::array<std::array<double,4>,3> matrix{};
                for(int i=0;i<3;++i) {
                    for(int j=0;j<3;++j)matrix[i][j]=value.jacobian[i][j];
                    matrix[i][3]=-value.residual[i];
                }
                for(int j=0;j<3;++j) {
                    int pivot=j;
                    for(int i=j+1;i<3;++i)if(std::abs(matrix[i][j])>std::abs(matrix[pivot][j]))pivot=i;
                    if(std::abs(matrix[pivot][j])<1.e-18)return false;
                    std::swap(matrix[j],matrix[pivot]);
                    const double diagonal=matrix[j][j];
                    for(int k=j;k<4;++k)matrix[j][k]/=diagonal;
                    for(int i=0;i<3;++i)if(i!=j) {
                        const double factor=matrix[i][j];
                        for(int k=j;k<4;++k)matrix[i][k]-=factor*matrix[j][k];
                    }
                }
                bool accepted=false;
                for(double damping=1.;damping>=1./65536.;damping*=.5) {
                    Reduced trial=z;for(int i=0;i<3;++i)trial[i]+=damping*matrix[i][3];
                    if(!std::isfinite(trial[0])||!std::isfinite(trial[1])||!std::isfinite(trial[2])||
                       trial[0]<=0.||trial[1]<=0.||trial[2]<0.)continue;
                    try {
                        // Eliminated linear pools need only be physical at the
                        // converged root; trial pressures and tissue O2 stay positive.
                        if(norm(reduced_evaluation(old,trial,dt,e,gas,false),trial)<n) {
                            z=trial;accepted=true;break;
                        }
                    } catch(const std::runtime_error&) {}
                }
                if(!accepted)return false;
            }
        } catch(const std::runtime_error&) {return false;}
        return false;
    }

    bool backward_euler(State& y,double dt,const Environment& e,const AlbaneseGas& gas) const {
        State trial=y;
        if(backward_euler_reduced(trial,dt,e,gas)) {y=trial;return true;}
        return backward_euler_full(y,dt,e,gas);
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
