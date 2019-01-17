#pragma once

#include "INIreader.h"
#include "gpumcd/Settings.h"

class DosiaSettings{
public:
	//from cmd args
	string rt_files;

	//from ini
	string gpumcd_material_data_dir;
	string gpumcd_machine_dir;
	string hounsfield_conversion_dir;

	int verbose;
	bool dbgoutput;
	
	float field_margin;
	bool dose_per_fraction;
	bool continous_materials;
	bool pinnacle_vmat_interpolation;
	bool monte_carlo_high_precision;
	bool score_dose_to_water;
	bool score_and_transport_in_water;
	bool in_aqua_vivo;
	
	bool gamma_comparison;
	bool gamma_global_dose;
	float gamma_isodose_region;
	float gamma_dd;
	float gamma_dta;

	string MRLinac_MV7;
	string Agility_MV6_FF;
	string Agility_MV6_NoFF;
	string Agility_MV10_FF;
	string Agility_MV10_NoFF;

	PhysicsSettings physicsSettings;
	PlanSettings planSettings;

	DosiaSettings() = default;
	DosiaSettings(const string &, const string &);
	DosiaSettings(const DosiaSettings &) = default;
	~DosiaSettings() = default;
};

DosiaSettings::DosiaSettings(const string &ini_file, const string & _rt_files) : rt_files(_rt_files){
	INIReader ini(ini_file);
	if (ini.ParseError() < 0) {
		throw std::pair<int, string>(1, "Can't load " + ini_file + "\n");
	}
	
	gpumcd_material_data_dir = ini.Get("directories", "gpumcd_material_data_dir", ".");
	gpumcd_machine_dir = ini.Get("directories", "gpumcd_machine_dir", ".");
	hounsfield_conversion_dir = ini.Get("directories", "hounsfield_conversion_dir", ".");

	MRLinac_MV7 = ini.Get("gpumcd_machines", "MRLinac_MV7", ".");
	Agility_MV6_FF = ini.Get("gpumcd_machines", "Agility_MV6_FF", ".");
	Agility_MV6_NoFF = ini.Get("gpumcd_machines", "Agility_MV6_NoFF", ".");
	Agility_MV10_FF = ini.Get("gpumcd_machines", "Agility_MV10_FF", ".");
	Agility_MV10_NoFF = ini.Get("gpumcd_machines", "Agility_MV10_NoFF", ".");

	verbose = ini.GetInteger("debug", "verbose", 0);
	dbgoutput = ini.GetBoolean("debug", "output", false);

	field_margin = ini.GetReal("dose", "field_margin", 5.f);
	dose_per_fraction = ini.GetBoolean("dose", "dose_per_fraction", true);
	continous_materials = ini.GetBoolean("dose", "continous_materials", true);
	pinnacle_vmat_interpolation = ini.GetBoolean("dose", "pinnacle_vmat_interpolation", false);
	monte_carlo_high_precision = ini.GetBoolean("dose", "monte_carlo_high_precision", false);
	score_dose_to_water = ini.GetBoolean("dose", "score_dose_to_water", false);
	score_and_transport_in_water = ini.GetBoolean("dose", "score_and_transport_in_water", false);
	in_aqua_vivo = ini.GetBoolean("dose", "in_aqua_vivo", false);

	gamma_comparison = ini.GetBoolean("gamma", "comparison", false);
	gamma_global_dose = ini.GetBoolean("gamma", "global_dose", true);
	gamma_isodose_region = ini.GetReal("gamma", "isodose_region", 10);
	gamma_dd = ini.GetReal("gamma", "dd", 3);
	gamma_dta = ini.GetReal("gamma", "dta", 2);
	
	// Load settings gpumcd physics, defaults from gpumcd/Settings.h
	physicsSettings.photonTransportCutoff = ini.GetReal("gpumcd_physicssettings", "photonTransportCutoff", 0.01f);
	physicsSettings.electronTransportCutoff = ini.GetReal("gpumcd_physicssettings", "electronTransportCutoff", 0.189f);
	physicsSettings.inputMaxStepLength = ini.GetReal("gpumcd_physicssettings", "inputMaxStepLength", 0.75f);
	physicsSettings.referenceMedium = ini.GetInteger("gpumcd_physicssettings", "referenceMedium", -1);
	physicsSettings.useElectronInAirSpeedup = ini.GetBoolean("gpumcd_physicssettings", "useElectronInAirSpeedup", true);
	physicsSettings.electronInAirSpeedupDensityThreshold = ini.GetReal("gpumcd_physicssettings", "electronInAirSpeedupDensityThreshold", 0.002f);

	// Load settings gpumcd plan, NOT ALL defaults from gpumcd/Settings.h, parts may be overridden by the RTBeam
	planSettings.goalSfom = ini.GetReal("gpumcd_plansettings", "goalSfom", 54.32f); //54.32f is a flag value, will be used to identify UNSET
	planSettings.statThreshold = ini.GetReal("gpumcd_plansettings", "statThreshold", 0.5f);
	planSettings.maxNumParticles = static_cast<uint64_t>(ini.GetReal("gpumcd_plansettings", "maxNumParticles", 1e10)); //GetReal so that we can cope with scientific notation
	planSettings.densityThresholdSfom = ini.GetReal("gpumcd_plansettings", "densityThresholdSfom", 0.2f);
	planSettings.densityThresholdOutput = ini.GetReal("gpumcd_plansettings", "densityThresholdOutput", 0.f);
	planSettings.useApproximateStatistics = ini.GetBoolean("gpumcd_plansettings", "useApproximateStatistics", true); // is AUTOMATICALLY set to FALSE for SEGMENT calcs.

	// Done. Print settings for log?
	if (verbose > 0){
		cerr << "Verbose level set at " << verbose << ".\n";

		cerr << "field_margin = " << field_margin << ".\n";
		cerr << "dose_per_fraction = " << dose_per_fraction << ".\n";
		cerr << "pinnacle_vmat_interpolation = " << pinnacle_vmat_interpolation << ".\n";
		cerr << "monte_carlo_high_precision = " << monte_carlo_high_precision << ".\n";

		if (gamma_comparison) cerr << "Gamma comparison enabled.\n";
		if (dbgoutput) cerr << "Debug outputs will be written to disk.\n";
		//if (in_aqua_vivo) cerr << "Forcing all densities inside patient threshold to 1.0g/cm3 (as EpidTrial.py in Pinnacle).\n";
		if (in_aqua_vivo) cerr << "in_aqua_vivo currently not correctly implemented. Will be removed. Dosia dump should fix this.\n";
		if (in_aqua_vivo) in_aqua_vivo = false;

		if (score_and_transport_in_water) { cerr << "Computing dose and transport in water instead of medium.\n"; }
		else if (score_dose_to_water) { cerr << "Computing dose to water instead of medium.\n"; };
	}
};