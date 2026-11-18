#ifndef ALBANESE_GAS_H
#define ALBANESE_GAS_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// Albanese et al. (2016), doi:10.1152/ajpheart.00230.2014, A46-A49.
// This supplies local equilibrium targets, not new transported quantities.
struct AlbaneseGas {
    struct Curve {
        double saturation_mmol_l, h, alpha, beta, K_mmHg;

        struct Gradient { double pressure, other; };
        Gradient gradient(double value, double pressure_mmHg, double other_mmHg,
                          double molar_volume_l_mol) const {
            const double maximum=saturation_mmol_l*1.e-3*molar_volume_l_mol;
            const double common=value*(1.-value/maximum)/h;
            return {common/pressure_mmHg,
                    common*(beta/(1.+beta*other_mmHg)-alpha/(1.+alpha*other_mmHg))};
        }

        double total(double pressure_mmHg, double other_mmHg,
                     double molar_volume_l_mol) const {
            if (!std::isfinite(pressure_mmHg) || pressure_mmHg < 0. ||
                !std::isfinite(other_mmHg) || other_mmHg < 0.)
                throw std::runtime_error("Albanese: invalid partial pressure");
            if (pressure_mmHg == 0.) return 0.;
            const double log_x = std::log(pressure_mmHg) - std::log(K_mmHg)
                + std::log1p(beta * other_mmHg) - std::log1p(alpha * other_mmHg);
            const double z = log_x / h; // The exponent in A46/A48 is 1/h.
            const double fraction = z >= 0. ? 1. / (1. + std::exp(-z))
                : std::exp(z) / (1. + std::exp(z));
            // mmol/l blood -> l(STPD)/l blood = m3(STPD)/m3 blood.
            return saturation_mmol_l * 1.e-3 * molar_volume_l_mol * fraction;
        }
    };

    bool enabled = false;
    Curve oxygen{9., .3836, .03198, .008275, 14.99};
    Curve carbon{86.11, 1.819, .05591, .03255, 194.4};
    double hb_co2_fraction = .05; // Additional component-partition assumption.
    double molar_volume_l_mol = 22.414; // STPD, not BTPS.
    double reference_hb_capacity = .201; // 15 g/dl * 1.34 ml O2/g / 100.

    static double number(const std::string& text) {
        std::size_t end = 0;
        const double value = std::stod(text, &end);
        if (end != text.size() || !std::isfinite(value))
            throw std::runtime_error("Albanese: invalid numeric parameter: " + text);
        return value;
    }

    bool read_row(const std::vector<std::string>& row) {
        if (row.empty()) return false;
        if (row[0] == "gas_model") {
            if (row.size() != 2 || (row[1] != "legacy" && row[1] != "albanese2016"))
                throw std::runtime_error("Expected gas_model,legacy|albanese2016");
            enabled = row[1] == "albanese2016";
        } else if (row[0] == "albanese_o2" || row[0] == "albanese_co2") {
            if (row.size() != 6)
                throw std::runtime_error("Expected albanese_o2/co2,Csat_mmol_l,h,alpha,beta,K_mmHg");
            Curve& curve = row[0] == "albanese_o2" ? oxygen : carbon;
            curve = {number(row[1]), number(row[2]), number(row[3]), number(row[4]), number(row[5])};
        } else if (row[0] == "albanese_options") {
            if (row.size() != 4)
                throw std::runtime_error("Expected albanese_options,f_Hb,molar_volume_l_mol,reference_hb_capacity");
            hb_co2_fraction = number(row[1]);
            molar_volume_l_mol = number(row[2]);
            reference_hb_capacity = number(row[3]);
        } else return false;
        validate();
        return true;
    }

    void validate() const {
        for (const Curve* c : {&oxygen, &carbon}) {
            for (double value : {c->saturation_mmol_l, c->h, c->alpha, c->beta, c->K_mmHg})
                if (!std::isfinite(value)) throw std::runtime_error("Albanese: nonfinite curve parameter");
            if (c->saturation_mmol_l <= 0. || c->h <= 0. || c->K_mmHg <= 0. ||
                c->alpha < 0. || c->beta < 0.)
                throw std::runtime_error("Albanese: invalid dissociation curve parameter");
        }
        if (!std::isfinite(hb_co2_fraction) || hb_co2_fraction < 0. || hb_co2_fraction > 1. ||
            !std::isfinite(molar_volume_l_mol) || molar_volume_l_mol <= 0. ||
            !std::isfinite(reference_hb_capacity) || reference_hb_capacity <= 0.)
            throw std::runtime_error("Albanese: invalid partition/volume/capacity parameter");
    }

    struct OxygenGradient { double pressure=0., other=0., plasma=0.; };
    double saturation(double plasma_o2, double po2_mmHg, double pco2_mmHg,
                      double capacity, OxygenGradient* derivative=nullptr) const {
        if(derivative) *derivative=OxygenGradient{};
        if (!std::isfinite(capacity) || capacity < 0. ||
            !std::isfinite(plasma_o2) || plasma_o2 < 0.)
            throw std::runtime_error("Albanese: invalid oxygen concentration/capacity");
        if (capacity == 0.) return 0.;
        // Match the empirical reference-Hb curve to the existing RBC*Z capacity.
        const double reference = oxygen.total(po2_mmHg, pco2_mmHg, molar_volume_l_mol);
        const double total = reference * capacity / reference_hb_capacity;
        const double raw=(total-plasma_o2)/capacity;
        if(derivative && raw>0. && raw<1.) {
            const auto g=oxygen.gradient(reference,po2_mmHg,pco2_mmHg,molar_volume_l_mol);
            *derivative={g.pressure/reference_hb_capacity,g.other/reference_hb_capacity,-1./capacity};
        }
        // Bound the saturation when extrapolating the empirical curve.
        return std::max(0., std::min(1., raw));
    }

    struct CarbonTargets { double bicarbonate, carbamino; };
    struct CarbonGradient { double bicarbonate_po, bicarbonate_pc, carbamino_po, carbamino_pc; };
    CarbonTargets carbon_targets(double po2_mmHg, double pco2_pa,
                                 double alpha_pla, double alpha_rbc,
                                 double phi_pla, double phi_rbc,
                                 double pa_per_mmHg, CarbonGradient* derivative=nullptr) const {
        const double total = carbon.total(pco2_pa / pa_per_mmHg, po2_mmHg, molar_volume_l_mol);
        const double dissolved = (phi_pla * alpha_pla + phi_rbc * alpha_rbc) * pco2_pa;
        const double bicarbonate = (1. - hb_co2_fraction) * total - dissolved;
        if (!std::isfinite(bicarbonate) || bicarbonate < -1.e-12)
            throw std::runtime_error("Albanese: negative equilibrium bicarbonate; check parameters/pressure range");
        if(derivative) {
            const auto g=carbon.gradient(total,pco2_pa/pa_per_mmHg,po2_mmHg,molar_volume_l_mol);
            *derivative={bicarbonate>0.?(1.-hb_co2_fraction)*g.other:0.,
                         bicarbonate>0.?(1.-hb_co2_fraction)*g.pressure/pa_per_mmHg
                            -(phi_pla*alpha_pla+phi_rbc*alpha_rbc):0.,
                         hb_co2_fraction*g.other,hb_co2_fraction*g.pressure/pa_per_mmHg};
        }
        return {std::max(0., bicarbonate), hb_co2_fraction * total};
    }
};

#endif
