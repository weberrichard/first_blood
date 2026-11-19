#include "../source/moc_node.h"
#include "../source/moc_edge.h"
#include <stdexcept>
#include <iostream>
#include <type_traits>
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
int main(){
    static_assert(std::is_same<decltype(std::declval<const moc_edge&>().get_velocity()),const std::vector<double>&>::value,"Velocity getter must not copy");
    static_assert(std::is_same<decltype(std::declval<const moc_edge&>().get_area()),const std::vector<double>&>::value,"Area getter must not copy");
    static_assert(std::is_same<decltype(std::declval<const moc_edge&>().get_new_velocity()),const std::vector<double>&>::value,"New velocity getter must not copy");
    moc_node node("test");node.type="junction";node.resistance=1.e10;
    node.is_resistance=1.;node.pressure_out=1.e5;node.density=1060.;
    node.do_save_memory=false;
    node.initialization(101000.,4.9e15,.9,.003,.026,.02,.7,.3,.05);
    node.boundary_variables(102000.,.001);
    require(node.current_pressure==102000.,"Current pressure was tied to saving");
    require(node.pressure.size()==1 && node.pressure.back()==101000.,"Disabled history grew");
    require(std::abs(node.current_flow_rate-2.e-7)<1.e-20,"Current leak flow is wrong");
    node.do_save_memory=true;node.boundary_variables(103000.,.002);
    require(node.pressure.back()==node.current_pressure && node.pressure.size()==2,"Enabled history differs from current state");
    node.save_field_variables(.003,104000.,.0001);
    require(node.current_pressure==104000. && node.current_flow_rate==.0001,"Prescribed current state was not updated");
    node.initialization(101000.,4.9e15,.9,.003,.026,.02,.7,.3,.05);
    require(node.current_pressure==101000. && node.current_flow_rate==0.,"Reinitialization retained stale state");
    double worst_flow=0.,worst_derivative=0.;
    for(int material:{0,1}) {
        moc_edge edge("test");edge.name="test";edge.length=.1;
        edge.nominal_diameter_start=.01;edge.nominal_diameter_end=.009;
        edge.nominal_thickness_start=.001;edge.nominal_thickness_end=.0009;
        edge.elasticity=4.e5;edge.resistance_start=0.;edge.resistance_end=0.;
        edge.division_points=8;edge.gravity=9.81;edge.density=1060.;
        edge.kinematic_viscosity=3.5e-6;edge.atmospheric_pressure=1.e5;
        edge.poisson_coefficient=.5;edge.courant_number=.5;
        edge.material_const={2.e6,-2250.,86500.};edge.do_save_memory=false;
        edge.initialization(1.e5,material,4.9e15,.9,.003,.026,.02,.7,.3,.05);
        for(bool start:{false,true}) {
            auto context=edge.prepare_junction(start,1.e-5);
            for(int i=-10;i<=10;++i) {
                const double pressure=1.e5+i*300.;
                const auto fast=context.evaluate(pressure);
                const auto original=start?edge.boundary_junction_start(pressure,1.e-5):edge.boundary_junction_end(pressure,1.e-5);
                const double h=.05;
                const double plus=(start?edge.boundary_junction_start(pressure+h,1.e-5):edge.boundary_junction_end(pressure+h,1.e-5))[0];
                const double minus=(start?edge.boundary_junction_start(pressure-h,1.e-5):edge.boundary_junction_end(pressure-h,1.e-5))[0];
                const double numeric=(plus-minus)/(2.*h);
                worst_flow=std::max(worst_flow,std::abs(fast[0]-original[0]));
                worst_derivative=std::max(worst_derivative,std::abs(fast[1]-numeric)/std::max(1.e-15,std::abs(numeric)));
            }
        }
    }
    require(worst_flow<1.e-13,"Prepared boundary changes original flow law");
    require(worst_derivative<1.e-6,"Analytic boundary derivative differs from central differences");
    std::cout<<"Hydraulic optimization: state/history and linear/Olufsen boundary checks passed; max flow error="<<worst_flow<<", derivative relative error="<<worst_derivative<<"\n";
}
