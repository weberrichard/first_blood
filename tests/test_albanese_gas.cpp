#include "../source/albanese_gas.h"
#include <iostream>

static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

int main() {
    AlbaneseGas gas;
    require(!gas.enabled, "Legacy must remain the default");
    gas.read_row({"gas_model", "albanese2016"});
    require(gas.enabled, "Model selector not applied");
    // At X=1 the content must be half saturation, with mmol/l converted to v/v.
    require(std::abs(gas.oxygen.total(14.99, 0., 22.414) - .100863) < 1.e-12,
            "Oxygen units/reference point incorrect");
    const double oxygen_low_co2 = gas.oxygen.total(60., 30., gas.molar_volume_l_mol);
    const double oxygen_high_co2 = gas.oxygen.total(60., 60., gas.molar_volume_l_mol);
    require(oxygen_low_co2 > oxygen_high_co2, "Bohr response has wrong direction");
    require(gas.carbon.total(40., 30., gas.molar_volume_l_mol) >
            gas.carbon.total(40., 100., gas.molar_volume_l_mol),
            "Haldane response has wrong direction");
    const double capacity = 4.9e15 * 3.73e-17;
    const double sat = gas.saturation(100.*3.11e-5, 100., 40., capacity);
    require(sat > .9 && sat < 1., "Reference Hb capacity is incompatible");
    const double d_hb = .9-sat;
    const double d_plasma = capacity*d_hb;
    require(std::abs(d_plasma/.4+capacity*(-d_hb/.4))<1.e-14,
            "Internal oxygen source creates gas");
    auto target = gas.carbon_targets(100., 40.*133.3616, 4.5e-6, 5.e-6,
                                    .56, .44, 133.3616);
    const double dissolved=(.56*4.5e-6+.44*5.e-6)*40.*133.3616;
    require(std::abs(target.bicarbonate+target.carbamino+dissolved-
                    gas.carbon.total(40.,100.,22.414))<1.e-12,
            "Carbon partition does not reconstruct total content");
    const double eta_b=target.bicarbonate-.44*.09-.56*.7;
    const double eta_h=target.carbamino-.44*.04;
    require(std::abs(.44*(-eta_b/.7-eta_h/.9)+.44*eta_b/.7+.44*eta_h/.9)<1.e-14,
            "Internal carbon reactions create gas");
    for (auto row : std::vector<std::vector<std::string>>{
             {"gas_model","unknown"}, {"albanese_o2","9","0","0","0","15"},
             {"albanese_options","1.1","22.414",".201"},
             {"albanese_co2","86","1.8","nan","0","194"},
             {"albanese_o2","9"}}) {
        bool rejected=false;
        try { AlbaneseGas invalid; invalid.read_row(row); }
        catch (const std::exception&) { rejected=true; }
        require(rejected,"Invalid parameter was accepted");
    }
    gas.read_row({"albanese_options",".08","22.4",".2"});
    require(gas.hb_co2_fraction==.08 && gas.molar_volume_l_mol==22.4,
            "File parameters were not applied");
    std::cout << "Albanese curve, coupling, units, capacity, conservation and input checks passed\n";
}
