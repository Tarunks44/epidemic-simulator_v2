//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#include "wastewater.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cassert>

using std::vector;
using std::string;

// ============================================================================
// FEATURE 1: RNA VIRAL LOAD CALCULATION
// ============================================================================

/**
 * Calculate viral shedding for an individual based on days since infection
 * Uses a log-normal-like curve peaking around day 5-7
 *
 * Shedding profile:
 * - Starts during exposed phase (day 1-2)
 * - Peaks during pre-symptomatic/early symptomatic (day 5-7)
 * - Continues into recovery (up to 20-30 days)
 * - Asymptomatic cases shed at ~60% of symptomatic levels
 */
double calculate_viral_shedding(double days_since_infection,
                               bool is_symptomatic,
                               const viral_shedding_params& params) {
    // No shedding before infection starts
    if (days_since_infection < params.shedding_start_day) {
        return 0.0;
    }

    // No shedding after duration ends
    if (days_since_infection > params.shedding_duration_days) {
        return 0.0;
    }

    // Log-normal-like shedding curve
    // Peak at peak_shedding_day, controlled by sigma
    double exponent = -0.5 * std::pow((days_since_infection - params.peak_shedding_day) /
                                      params.shedding_curve_sigma, 2);
    double shedding = params.peak_shedding_rate * std::exp(exponent);

    // Asymptomatic individuals shed less
    if (!is_symptomatic) {
        shedding *= params.asymptomatic_shedding_fraction;
    }

    return shedding;
}

/**
 * Add measurement noise to simulate real wastewater sampling variability
 * Typical CV (coefficient of variation) for WW measurements: 10-20%
 */
double add_measurement_noise(double concentration, double noise_std_fraction) {
    if (concentration <= 0) {
        return 0.0;
    }

    // Add Gaussian noise with std = noise_std_fraction * concentration
    double noise = uniform_real(-noise_std_fraction, noise_std_fraction);
    double noisy_concentration = concentration * (1.0 + noise);

    // Ensure non-negative
    return std::max(0.0, noisy_concentration);
}

/**
 * Check if a sample should be collected today based on sampling frequency
 */
bool should_collect_sample(const sewershed& ss, int current_day) {
    if (!ss.sampling_active) {
        return false;
    }

    // Calculate sampling interval in days
    double days_per_sample = 7.0 / ss.sampling_frequency_per_week;

    // Check if enough days have passed since last sample
    int days_since_last_sample = current_day - ss.last_sample_day;

    return (days_since_last_sample >= static_cast<int>(days_per_sample));
}

// ============================================================================
// INITIALIZATION
// ============================================================================

/**
 * Initialize sewersheds by mapping homes to sewershed catchment areas
 *
 * Default strategy: Group wards into sewersheds
 * - Each sewershed serves multiple wards
 * - Population and flow rates are calculated from homes
 */
void initialize_sewersheds(vector<sewershed>& sewersheds,
                          const vector<house>& homes,
                          const vector<community>& communities,
                          const wastewater_config& config) {

    if (!config.wastewater_enabled) {
        return;
    }

    // Clear existing sewersheds
    sewersheds.clear();

    // For now, create sewersheds based on ward groupings
    // You can customize this based on your city's sewage infrastructure

    // Group wards into sewersheds (example: 5-10 wards per sewershed)
    int wards_per_sewershed = 10;
    count_type total_wards = GLOBAL.num_wards;
    int num_sewersheds = (total_wards + wards_per_sewershed - 1) / wards_per_sewershed;

    for (int ss_id = 0; ss_id < num_sewersheds; ++ss_id) {
        sewershed ss;
        ss.sewershed_id = ss_id;
        ss.sewershed_name = "Sewershed_" + std::to_string(ss_id);
        ss.sampling_active = true;
        ss.sampling_frequency_per_week = 3.0; // Default: 3 samples per week
        ss.collection_efficiency = 0.85;
        ss.dilution_factor = 1.0;
        ss.last_sample_day = -10; // Allow immediate sampling

        // Assign wards to this sewershed
        int start_ward = ss_id * wards_per_sewershed;
        int end_ward = std::min(start_ward + wards_per_sewershed, static_cast<int>(total_wards));

        for (int ward = start_ward; ward < end_ward; ++ward) {
            ss.ward_ids.push_back(ward);
        }

        sewersheds.push_back(ss);
    }

    // Map homes to sewersheds and calculate population/flow rates
    for (count_type home_idx = 0; home_idx < homes.size(); ++home_idx) {
        const house& h = homes[home_idx];
        int ward_id = h.community; // Ward index

        // Find which sewershed this ward belongs to
        for (auto& ss : sewersheds) {
            if (std::find(ss.ward_ids.begin(), ss.ward_ids.end(), ward_id) != ss.ward_ids.end()) {
                ss.home_indices.push_back(home_idx);
                ss.population += h.individuals.size();
                break;
            }
        }
    }

    // Calculate flow rates based on population
    // Typical: 150-200 L/person/day for domestic wastewater
    double liters_per_person_per_day = 175.0;

    for (auto& ss : sewersheds) {
        ss.flow_rate_L_per_day = ss.population * liters_per_person_per_day;

        std::cout << "Initialized " << ss.sewershed_name
                  << ": Population=" << ss.population
                  << ", Wards=" << ss.ward_ids.size()
                  << ", Homes=" << ss.home_indices.size()
                  << ", Flow=" << ss.flow_rate_L_per_day << " L/day\n";
    }
}

// ============================================================================
// MAIN UPDATE FUNCTION
// ============================================================================

/**
 * Update wastewater surveillance - called each simulation step
 * Implements all three features:
 * 1. Calculate RNA viral load
 * 2. Aggregate detected clinical cases
 * 3. Estimate detection ratio
 */
void update_wastewater_surveillance(vector<sewershed>& sewersheds,
                                   const vector<agent>& agents,
                                   const vector<house>& homes,
                                   int current_time,
                                   const wastewater_config& config) {

    if (!config.wastewater_enabled || sewersheds.empty()) {
        return;
    }

    int current_day = current_time / GLOBAL.SIM_STEPS_PER_DAY;

    // Reset daily counters for all sewersheds
    for (auto& ss : sewersheds) {
        ss.reset_daily_counters();
    }

    // ========================================================================
    // FEATURE 1: Calculate RNA viral load in wastewater
    // ========================================================================

    for (auto& ss : sewersheds) {
        double total_viral_shedding = 0.0;

        // Iterate through all homes in this sewershed
        for (int home_idx : ss.home_indices) {
            const house& h = homes[home_idx];

            // Iterate through all individuals in this home
            for (int agent_idx : h.individuals) {
                const agent& person = agents[agent_idx];

                // Check if person is infected (not susceptible, recovered, or dead)
                if (person.infection_status != Progression::susceptible &&
                    person.infection_status != Progression::recovered &&
                    person.infection_status != Progression::dead) {

                    // Calculate days since infection
                    double days_infected = current_time - person.time_of_infection;
                    days_infected = days_infected / GLOBAL.SIM_STEPS_PER_DAY;

                    // Calculate viral shedding for this person
                    bool is_symptomatic = person.entered_symptomatic_state;
                    double shedding = calculate_viral_shedding(days_infected,
                                                              is_symptomatic,
                                                              config.shedding_params);

                    total_viral_shedding += shedding;
                }
            }
        }

        // Store total viral load
        ss.daily_total_viral_shedding = total_viral_shedding;

        // Calculate concentration in wastewater (copies/L)
        if (ss.flow_rate_L_per_day > 0) {
            ss.daily_wastewater_concentration =
                (total_viral_shedding * ss.collection_efficiency) /
                (ss.flow_rate_L_per_day * ss.dilution_factor);
        }

        // Normalized concentration per 100,000 population
        if (ss.population > 0) {
            ss.normalized_concentration =
                (ss.daily_wastewater_concentration * 100000.0) / ss.population;
        }
    }

    // ========================================================================
    // FEATURE 2: Aggregate detected clinical cases by sewershed
    // ========================================================================

    aggregate_clinical_cases_by_sewershed(sewersheds, agents, homes, current_time);

    // ========================================================================
    // FEATURE 3: Estimate detection ratio
    // ========================================================================

    estimate_detection_ratio(sewersheds, config);
}

// ============================================================================
// FEATURE 2: AGGREGATE CLINICAL CASES BY SEWERSHED
// ============================================================================

/**
 * Aggregate clinical surveillance data (detected cases, hospitalizations, deaths)
 * by sewershed area
 */
void aggregate_clinical_cases_by_sewershed(vector<sewershed>& sewersheds,
                                          const vector<agent>& agents,
                                          const vector<house>& homes,
                                          int current_time) {

    int current_day = current_time / GLOBAL.SIM_STEPS_PER_DAY;

    for (auto& ss : sewersheds) {
        count_type detected_today = 0;
        count_type hospitalised_today = 0;
        count_type deaths_today = 0;

        // Iterate through all homes in this sewershed
        for (int home_idx : ss.home_indices) {
            const house& h = homes[home_idx];

            // Iterate through all individuals in this home
            for (int agent_idx : h.individuals) {
                const agent& person = agents[agent_idx];

                // Count detected cases (using case_infection_ratio)
                // A case is "detected" when they become symptomatic
                if (person.entered_symptomatic_state) {
                    // Check if they became symptomatic today
                    double days_since_infection = current_time - person.time_of_infection;
                    double days_to_symptomatic = (person.incubation_period + person.asymptomatic_period);

                    // If they just became symptomatic (within last timestep)
                    if (std::abs(days_since_infection - days_to_symptomatic) < 1.5) {
                        // Apply detection probability
                        if (bernoulli(GLOBAL.case_infection_ratio)) {
                            detected_today++;
                        }
                    }
                }

                // Count hospitalizations
                if (person.entered_hospitalised_state) {
                    double days_since_infection = current_time - person.time_of_infection;
                    double days_to_hospitalised = (person.incubation_period +
                                                   person.asymptomatic_period +
                                                   person.symptomatic_period);

                    if (std::abs(days_since_infection - days_to_hospitalised) < 1.5) {
                        hospitalised_today++;
                    }
                }

                // Count deaths
                if (person.infection_status == Progression::dead) {
                    // This is simplistic - you may want to track when death occurred
                    deaths_today++;
                }
            }
        }

        // Update sewershed counters
        ss.daily_detected_cases = detected_today;
        ss.daily_hospitalised = hospitalised_today;
        ss.daily_deaths = deaths_today;

        // Update cumulative counters
        ss.cumulative_detected_cases += detected_today;
        ss.cumulative_hospitalised += hospitalised_today;
        ss.cumulative_deaths += deaths_today;
    }
}

// ============================================================================
// FEATURE 3: DETECTION RATIO ESTIMATION
// ============================================================================

/**
 * Estimate true infections from wastewater signal and calculate detection ratio
 *
 * Detection ratio = detected_cases / estimated_true_infections
 *
 * This tells us what fraction of infections are being detected by clinical surveillance
 */
void estimate_detection_ratio(vector<sewershed>& sewersheds,
                             const wastewater_config& config) {

    for (auto& ss : sewersheds) {
        // Estimate true infections from wastewater concentration
        // This is a simplified model - in practice, you'd calibrate this

        // Method: Higher viral concentration -> more infections
        // Calibration factor converts copies/L to estimated infections per 100k

        if (ss.daily_wastewater_concentration > config.shedding_params.detection_limit_copies_per_L) {
            // Estimate infections per 100k population
            double infections_per_100k = ss.daily_wastewater_concentration *
                                        config.viral_load_to_infection_factor;

            // Scale to actual population
            ss.estimated_true_infections = static_cast<count_type>(
                (infections_per_100k * ss.population) / 100000.0
            );
        } else {
            // Below detection limit
            ss.estimated_true_infections = 0;
        }

        // Calculate detection ratio
        // detection_ratio = detected / true
        if (ss.estimated_true_infections > config.minimum_infections_for_detection) {
            if (ss.cumulative_detected_cases > 0) {
                ss.detection_ratio = static_cast<double>(ss.cumulative_detected_cases) /
                                    static_cast<double>(ss.estimated_true_infections);
            } else {
                ss.detection_ratio = 0.0;
            }

            // Calculate underreporting factor (how many true cases per detected case)
            if (ss.detection_ratio > 0.01) {
                ss.underreporting_factor = 1.0 / ss.detection_ratio;
            } else {
                ss.underreporting_factor = 0.0;
            }
        } else {
            ss.detection_ratio = 0.0;
            ss.underreporting_factor = 0.0;
        }
    }
}

// ============================================================================
// GROUND TRUTH CALCULATION (for validation)
// ============================================================================

/**
 * Get the true number of active infections in a sewershed
 * This is the ground truth from the simulation (not available in real world)
 */
count_type get_true_active_infections_in_sewershed(const sewershed& ss,
                                                   const vector<agent>& agents,
                                                   const vector<house>& homes) {
    count_type true_infections = 0;

    for (int home_idx : ss.home_indices) {
        const house& h = homes[home_idx];

        for (int agent_idx : h.individuals) {
            const agent& person = agents[agent_idx];

            // Count all infected (exposed, infective, symptomatic, hospitalised, critical)
            if (person.infection_status != Progression::susceptible &&
                person.infection_status != Progression::recovered &&
                person.infection_status != Progression::dead) {
                true_infections++;
            }
        }
    }

    return true_infections;
}

// ============================================================================
// SAMPLE GENERATION
// ============================================================================

/**
 * Generate a complete wastewater sample record for output
 */
wastewater_sample generate_sample_record(const sewershed& ss,
                                        int current_day,
                                        const vector<agent>& agents,
                                        const vector<house>& homes) {
    wastewater_sample sample;

    sample.sewershed_id = ss.sewershed_id;
    sample.day = current_day;

    // Feature 1: RNA load
    sample.viral_concentration_copies_per_L = ss.daily_wastewater_concentration;
    sample.normalized_concentration_per_100k = ss.normalized_concentration;
    sample.total_viral_load_copies = ss.daily_total_viral_shedding;

    // Feature 2: Detected cases
    sample.detected_cases_in_area = ss.daily_detected_cases;
    sample.hospitalised_in_area = ss.daily_hospitalised;
    sample.deaths_in_area = ss.daily_deaths;
    sample.cumulative_detected = ss.cumulative_detected_cases;
    sample.cumulative_hospitalised = ss.cumulative_hospitalised;
    sample.cumulative_deaths = ss.cumulative_deaths;

    // Feature 3: Detection ratio
    sample.estimated_true_infections = ss.estimated_true_infections;
    sample.true_active_infections = get_true_active_infections_in_sewershed(ss, agents, homes);
    sample.detection_ratio = ss.detection_ratio;
    sample.underreporting_factor = ss.underreporting_factor;

    // Metadata
    sample.population = ss.population;
    sample.below_detection_limit = (ss.daily_wastewater_concentration < 1000.0);

    return sample;
}

// ============================================================================
// OUTPUT FUNCTIONS
// ============================================================================

/**
 * Output wastewater samples to CSV file
 */
void output_wastewater_samples(const string& output_path,
                              const vector<wastewater_sample>& samples) {

    std::ofstream fout(output_path, std::ios::out);

    if (!fout) {
        std::cerr << "Error: Could not open wastewater output file: " << output_path << "\n";
        return;
    }

    // Write CSV header
    fout << "Day,Sewershed_ID,Population,"
         << "Viral_Concentration_copies_per_L,Normalized_Conc_per_100k,Total_Viral_Load,"
         << "Detected_Cases,Hospitalised,Deaths,"
         << "Cumulative_Detected,Cumulative_Hospitalised,Cumulative_Deaths,"
         << "Estimated_True_Infections,True_Active_Infections,"
         << "Detection_Ratio,Underreporting_Factor,Below_Detection_Limit\n";

    // Write data
    for (const auto& sample : samples) {
        fout << sample.day << ","
             << sample.sewershed_id << ","
             << sample.population << ","
             << sample.viral_concentration_copies_per_L << ","
             << sample.normalized_concentration_per_100k << ","
             << sample.total_viral_load_copies << ","
             << sample.detected_cases_in_area << ","
             << sample.hospitalised_in_area << ","
             << sample.deaths_in_area << ","
             << sample.cumulative_detected << ","
             << sample.cumulative_hospitalised << ","
             << sample.cumulative_deaths << ","
             << sample.estimated_true_infections << ","
             << sample.true_active_infections << ","
             << sample.detection_ratio << ","
             << sample.underreporting_factor << ","
             << (sample.below_detection_limit ? 1 : 0) << "\n";
    }

    fout.close();
    std::cout << "Wrote " << samples.size() << " wastewater samples to " << output_path << "\n";
}

/**
 * Output sewershed timeseries data
 */
void output_sewershed_timeseries(const string& output_path,
                                const vector<vector<wastewater_sample>>& timeseries) {

    // timeseries[sewershed_id][time_idx]

    for (size_t ss_id = 0; ss_id < timeseries.size(); ++ss_id) {
        string filename = output_path + "/sewershed_" + std::to_string(ss_id) + "_timeseries.csv";
        output_wastewater_samples(filename, timeseries[ss_id]);
    }
}
