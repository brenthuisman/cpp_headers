#pragma once

//#include <set>
//#include <string>
//using std::string;
#include <assert.h>

#include "CalculationInformation.h"

#include "pystring.h"
using namespace pystring;
using std::string;
#include "tools.h"
using namespace parse;
using namespace vect;
#include "settings.h"

class RTBeam;
class Parser;
class dicomParser;
class pinnacleParser;

struct ControlPoint {
	// gpumcd computes one ModifierInformation+BeamInformation pair at a time, and a full description of a controlpoint consists of both
	ModifierInformation collimator; //ASMX == parallelJaw
	BeamInformation beamInfo;
};

enum class AcceleratorType { UNKNOWN, EMPTY, MLCi80, Agility, MRLinac };
enum class Energy { UNKNOWN, MV6, MV7, MV10 };
enum class Filter { UNKNOWN, FF, NoFF }; //YES is the normal condition.

class Accelerator {
public:
	AcceleratorType type;
	int leafs_per_bank;
	Energy energy;
	Filter filter = Filter::FF;//default is WITH flattening filter, is overridden to FFF if encountered

	Accelerator(AcceleratorType _type){
		if (_type == AcceleratorType::MLCi80){
			leafs_per_bank = 40;
		}
		else if (_type == AcceleratorType::Agility){
			leafs_per_bank = 80;
		}
		else if (_type == AcceleratorType::MRLinac){
			leafs_per_bank = 80;
		}
		type = _type;
	}
	Accelerator(){ type = AcceleratorType::EMPTY; };
};

enum class BeamType { UNKNOWN, IMRT, VMAT };

struct BeamMetaData {
	//plan variables
	BeamType beamtype;
	float weight;
	float mu_per_fraction;// = 1.f;
	int nr_fractions;// = 1;
	Accelerator accelerator;
	float prescriptiondose;
	string isocentername;
	Float3 isoc;
	Float3 bfield;

	//ct metadata
	float hu_slope;
	float hu_intercept;
	float outsidepatientairthreshold;
	int outsidepatientisctnumber;
	string patient_position;
	float couchremovalycoordinate;
	float couch_height;
	
	//params
	bool dose_per_fraction;
	float fieldMargin;
	bool pinnacle_vmat_interpolation;
	bool table_type;
};


class RTBeam{ //base
public:
	//members
	BeamMetaData metaData;
	vector<ControlPoint> controlPoints;

	//methods
	void printInfo();
	void printFirstLeaf();
	int num_cps(){ return controlPoints.size(); };

	//ctors
	RTBeam() = default;
	RTBeam(const RTBeam &) = default;
	~RTBeam() = default;

	RTBeam(DosiaSettings &);

private:
	DosiaSettings sett;
};


class Parser{
protected:
	//members
	BeamMetaData metaData;
	vector<ControlPoint> controlPoints;

	//params
	bool debug;
	
	//methods
	virtual void parseTrial(const string &){}; //default implementations do nothing
	virtual void parsePlan(const string &){};
	virtual void parseBeam(const string &){};
	virtual void parseScan(const string &){};
	virtual void parseDose(const string &){};
	virtual void setCPIs();
	virtual int num_cps(){ return controlPoints.size(); };

	//ctors, which cant be used by anyone
	Parser() = default;
	Parser(const Parser &) = default;
	~Parser() = default;
	Parser(BeamMetaData, bool = false); //metaData, debug (metaData moet pinnacleVMATmode, dose_per_fraction, fieldMargin geset hebben)

	friend class RTBeam; //only RTBeam can touch or construct this class.
};


class pinnacleParser : virtual private Parser{
public:

private:
	pinnacleParser(BeamMetaData, bool = false);
	void parseTrial(const string &);
	void parsePlan(const string &);
	void parseBeam(const string &);
	//void parseScan(const string &){};//not needed
	void parseDose(const string &);//prescription dose
	//void setCPIs();

	//ctors

	friend class RTBeam; //only RTBeam can touch or construct this class.
};


class dicomParser : virtual private Parser{
public:

private:
	dicomParser(BeamMetaData, bool = false);
	//void parseTrial(const string &){};//doesnt exist
	//void parsePlan(const string &); //TODO NOT IMPLEMENTED, SUSPICION: NO RELEVANT INFO HERE
	void parseBeam(const string &); //TODO
	void parseScan(const string &); //intercept, slope
	//void parseDose(const string &){};//prescription dose
	void setCPIs();

	/*
	'NumberOfControlPoints',   //BEAM_NRCONTROLPOINTS
		'BeamLimitingDeviceAngle', //BEAM_COLLIMATOR
		'GantryAngle',             //BEAM_GANTRY,
		'PatientSupportAngle',     //BEAM_COUCH
		'NominalBeamEnergy[0]'     //BEAM_ENERGY
		*/

	//ctors

	friend class RTBeam; //only RTBeam can touch or construct this class.
};


Parser::Parser(BeamMetaData _metaData, bool _debug) : debug(_debug), metaData(_metaData){
}


pinnacleParser::pinnacleParser(BeamMetaData _metaData, bool _debug) : Parser(_metaData, _debug){
}


dicomParser::dicomParser(BeamMetaData _metaData, bool _debug) : Parser(_metaData, _debug){
}


void pinnacleParser::parseTrial(const string &inFile) {
	auto dump = load_dump(inFile);
	assert(!dump.empty());

	for (auto &line : dump) {
		/*if (startswith(line.first, "requestedmonitorunitsperfraction[")) {	//better take this from dose.dump
			metaData.mu_per_fraction = stof(line.second);
			continue;
		}*/
		if (startswith(line.first, "numberoffractions[")) {
			metaData.nr_fractions = stoi(line.second);
			continue;
		}
		if (startswith(line.first, "negativemupenalty")) {
			metaData.hu_intercept = -stoi(line.second); //NEGATIVE!!!
			metaData.hu_slope = 1.f; //There is no slope in pinnacle.
			continue;
		}
		if (startswith(line.first, "outsidepatientairthreshold")) {
			metaData.outsidepatientairthreshold = stof(line.second);
			continue;
		}
		if (startswith(line.first, "outsidepatientisctnumber")) {
			metaData.outsidepatientisctnumber = stoi(line.second); // this is pre-offset
			continue;
		}
		if (startswith(line.first, "patient_position")) {
			metaData.patient_position = line.second;
			continue;
		}
		if (startswith(line.first, "couchremovalycoordinate")) {
			metaData.couchremovalycoordinate = stof(line.second);
			continue;
		}
	}

	if (debug) fprintf(stderr,"RTPLan: dicom intercept, slope: %.2f,%.2f\n", metaData.hu_intercept, metaData.hu_slope);
	if (debug) fprintf(stderr,"RTPLan: mu_per_fraction: %.2f.\nnr_fractions: %i\n", metaData.mu_per_fraction, metaData.nr_fractions);
}

void pinnacleParser::parseDose(const string &inFile) {
	auto dump = load_dump(inFile);
	assert(!dump.empty());

	for (auto &line : dump) {
		if (startswith(line.first, "prescriptiondose")) {
			metaData.prescriptiondose = stof(line.second);
			continue;
		}
		if (startswith(line.first, "requestedmonitorunitsperfraction")) {
			metaData.mu_per_fraction = stof(line.second);
			continue;
		}
	}
}

void dicomParser::parseScan(const string &inFile) {
	auto dump = load_dump(inFile);
	assert(!dump.empty());

	for (auto &line : dump) {
		if (startswith(line.first, "RescaleIntercept")) {
			metaData.hu_intercept = stof(line.second);
			continue;
		}
		if (startswith(line.first, "RescaleSlope")) {
			metaData.hu_slope = stof(line.second);
			continue;
		}
		if (startswith(line.first, "ImagePositionPatient")) {
			/*
			auto pos = split_vec<float>(line.second, "\\");
			fprintf(stderr, "RTPLan: dicom ImagePositionPatient : %.2f,%.2f,%.2f\n", x,y,z);
			*/
			continue;
		}

		
	}
	if (debug) fprintf(stderr, "RTPLan: dicom intercept, slope: %.2f,%.2f\n", metaData.hu_intercept, metaData.hu_slope);
}


void pinnacleParser::parsePlan(const string &inFile) {
	auto dump = load_dump(inFile);
	assert(!dump.empty());

	assert(!metaData.isocentername.empty());

	for (auto &line : dump) {
		if (startswith(lower(line.first), lower(metaData.isocentername) + "_x")) { //normalize case for isoc name
			metaData.isoc.x = stof(line.second);
			continue;
		}
		if (startswith(lower(line.first), lower(metaData.isocentername) + "_y")) { //minus!
			metaData.isoc.y = -stof(line.second);
			continue;
		}
		if (startswith(lower(line.first), lower(metaData.isocentername) + "_z")) { //minus!
			metaData.isoc.z = -stof(line.second);
			continue;
		}
	}
	if (debug) fprintf(stderr,"RTPLan: isoc %.2f %.2f %.2f\n", metaData.isoc.x, metaData.isoc.y, metaData.isoc.z);
}

void dicomParser::parseBeam(const string &inFile) {
	auto dump = load_dump(inFile);
	assert(!dump.empty());

	vector<int> ReferencedBeamNumber;
	vector<float> BeamDose;
	vector<float> BeamMeterset;

	vector<string> BeamLimitingDeviceSequence(3); //should usually have 2 or 3.

	for (auto &line : dump) {
		if (startswith(line.first, "FractionGroupSequence")) { //here we look for the beam weight, which is not yet known here
			vector<string> line_breakdown = split(line.first, ".");
			if (line_breakdown.size() <= 1) continue; //one or less: nothing to split. ignore
			
			if (line_breakdown[1] == "NumberOfFractionsPlanned"){
				metaData.nr_fractions = stoi(line.second);
			}

			int blds = get_index(line_breakdown[1]); //centerpiece has beam index
			if (blds < 0) continue; // can be no number, then skip

			if (line_breakdown[2] == "BeamDose"){
				BeamDose.push_back(stof(line.second));
			}
			if (line_breakdown[2] == "BeamMeterset"){
				BeamMeterset.push_back(stof(line.second));
			}
			continue;
		}
		if (line.first == "TreatmentMachineName") {
			if (startswith(line.second, "MLC160")){
				metaData.accelerator = Accelerator(AcceleratorType::Agility);
				continue;
			}
			else if (startswith(line.second, "M160")){
				metaData.accelerator = Accelerator(AcceleratorType::Agility);
				continue;
			}
			else if (startswith(line.second, "MLC80")){
				metaData.accelerator = Accelerator(AcceleratorType::MLCi80);
				continue;
			}
			else if (startswith(line.second, "M80")){
				metaData.accelerator = Accelerator(AcceleratorType::MLCi80);
				continue;
			}
			else {
				metaData.accelerator = Accelerator(AcceleratorType::UNKNOWN);
				continue;
			}
		}
		if (startswith(line.first, "BeamLimitingDeviceSequence")) {
			vector<string> line_breakdown = split(line.first, ".");
			if (line_breakdown.size() <= 1) continue; //one or less: nothing to split. ignore

			int blds = get_index(line_breakdown[0]);
			if (blds < 0) continue; // can be no number, then skip
			if (blds > 2) throw std::pair<int, string>(24, "Unexpected number of BeamLimitingDevices.");

			if (line_breakdown[1] == "RTBeamLimitingDeviceType"){
				BeamLimitingDeviceSequence[blds] = line.second;
			}
			if (line_breakdown[1] == "NumberOfLeafJawPairs"){
				if (BeamLimitingDeviceSequence[blds] == "MLCX"){
					if (metaData.accelerator.leafs_per_bank != stoi(line.second)) throw std::pair<int, string>(25, "Number of Leafs does not correspond to specified Accelerator.");
				}
			}
		}
		if (line.first == "BeamNumber") {
			assert(len(BeamDose) == len(BeamMeterset) == len(ReferencedBeamNumber));
			assert(len(ReferencedBeamNumber) > 0);
			//metaData.weight = BeamDose[stoi(line.first)]; // TODO: not sure which
			metaData.weight = BeamMeterset[index(ReferencedBeamNumber, stoi(line.first))];
			continue;
		}
		if (line.first == "NumberOfControlPoints") {
			controlPoints.resize(stoi(line.second));
			for (int i = 0; i < stoi(line.second); i++){
				controlPoints[i].collimator.mlc.leftLeaves.resize(metaData.accelerator.leafs_per_bank); // leaf_per_bank should be set before
				controlPoints[i].collimator.mlc.rightLeaves.resize(metaData.accelerator.leafs_per_bank);
			}
			continue;
		}
		if (startswith(line.first, "BeamType")) {
			if (startswith(line.second, "DYNAMIC")){
				metaData.beamtype = BeamType::VMAT;
			}
			//following copied from pinnacle::parseBeam
			else if (startswith(line.second, "Dynamic Arc")){
				metaData.beamtype = BeamType::VMAT;
			}
			else if (startswith(line.second, "Step & Shoot MLC")){
				metaData.beamtype = BeamType::IMRT;
			}
			else if (startswith(line.second, "Static")){
				metaData.beamtype = BeamType::IMRT;
			}
			else {
				metaData.beamtype = BeamType::UNKNOWN;
			}
			continue;
		}
		if (startswith(line.first, "PatientSetupSequence")) {
			vector<string> line_breakdown = split(line.first, ".");
			if (line_breakdown[1] == "PatientPosition"){
				metaData.patient_position = line.second;
			}
			continue;
		}
		if (startswith(line.first, "ControlPointSequence")) {
			vector<string> line_breakdown = split(line.first, ".");
			if (line_breakdown.size() <= 1) continue; //one or less: nothing to split. ignore

			//get CPI
			int cpi = get_index(line_breakdown[0]);
			if (cpi < 0) continue; // can be no number, then skip
			if (cpi > num_cps() - 1) throw std::pair<int, string>(26, "CPI out of range.");

			if (line_breakdown[1] == "ControlPointIndex") {
				assert(cpi == stoi(line.second));
				continue;
			}
			if (line_breakdown[1] == "NominalBeamEnergy") {
				if (startswith(line.second, "6")){
					metaData.accelerator.energy = Energy::MV6;
				}
				else if (startswith(line.second, "10")){
					metaData.accelerator.energy = Energy::MV10;
				}
				else if (startswith(line.second, "7")){ //cant happen I guess
					metaData.accelerator.energy = Energy::MV7;
				}
				else {
					metaData.accelerator.energy = Energy::UNKNOWN;
				}
				continue;
			}
			if (line_breakdown[1] == "NumberOfCompensators") { //TODO
				/*if (endswith(line.second, "FFF")){
				metaData.accelerator.filter = Filter::NoFF;
				//TODO: imrtfilter="Compensator" == FFF?
				}*/
				continue;
			}
			if (line_breakdown[1] == "GantryAngle") {
				controlPoints[cpi].beamInfo.gantryAngle = { stof(line.second), stof(line.second) };
				continue;
			}
			if (line_breakdown[1] == "PatientSupportAngle") {
				controlPoints[cpi].beamInfo.couchAngle = { stof(line.second), stof(line.second) };
				continue;
			}
			if (line_breakdown[1] == "BeamLimitingDeviceAngle") {
				controlPoints[cpi].beamInfo.collimatorAngle = { stof(line.second), stof(line.second) };
				continue;
			}
			if (line_breakdown[1] == "IsocenterPosition") {
				auto pos = types::split<float>(line.second, "\\");
				controlPoints[cpi].beamInfo.isoCenter.x = pos[0];
				controlPoints[cpi].beamInfo.isoCenter.y = pos[1];
				controlPoints[cpi].beamInfo.isoCenter.z = pos[2];
				continue;
			}
			if (line_breakdown[1] == "CumulativeMetersetWeight") {
				if (cpi == 0) {
					controlPoints[cpi].beamInfo.relativeWeight = stof(line.second);
				}
				else {
					controlPoints[cpi].beamInfo.relativeWeight = stof(line.second) - controlPoints[cpi - 1].beamInfo.relativeWeight;
				}
				continue;
			}
			if (startswith(line_breakdown[1], "BeamLimitingDevicePositionSequence")) {
				int bldi = get_index(line_breakdown[1]);
				if (bldi < 0) continue; // can be no number, then skip

				if (line_breakdown[2] == "LeafJawPositions") {
					auto pos = types::split<float>(line.second,"\\");
					if (BeamLimitingDeviceSequence[bldi] == "ASMX"){
						controlPoints[get_index(line.first)].collimator.parallelJaw.j1 = { pos[0], pos[0] }; //J1 always most negative coord according to doc.
						controlPoints[get_index(line.first)].beamInfo.fieldMin.first = pos[0] - metaData.fieldMargin;

						controlPoints[get_index(line.first)].collimator.parallelJaw.j2 = { pos[1], pos[1] };
						controlPoints[get_index(line.first)].beamInfo.fieldMax.first = pos[1] + metaData.fieldMargin;
					}
					if (BeamLimitingDeviceSequence[bldi] == "ASMY"){
						controlPoints[get_index(line.first)].collimator.perpendicularJaw.j1 = { pos[0], pos[0] };
						controlPoints[get_index(line.first)].beamInfo.fieldMin.second = pos[0] - metaData.fieldMargin;

						controlPoints[get_index(line.first)].collimator.perpendicularJaw.j2 = { pos[1], pos[1] };
						controlPoints[get_index(line.first)].beamInfo.fieldMax.second = pos[1] + metaData.fieldMargin;
					}
					if (BeamLimitingDeviceSequence[bldi] == "MLCX"){
						//pos heeft right bank eerst
						for (int p = 0; p < metaData.accelerator.leafs_per_bank; p++){
							controlPoints[cpi].collimator.mlc.rightLeaves[p] = { pos[p], pos[p] };
							controlPoints[cpi].collimator.mlc.leftLeaves[p] = { pos[p + metaData.accelerator.leafs_per_bank], pos[p + metaData.accelerator.leafs_per_bank] };
						}
						//NOOT: for square fields dicom does not necesarily give MLCX positions. Only jaw positions may be encountered!!
					}
				}
				continue;
			}
			continue;
		}
	}
}


void pinnacleParser::parseBeam(const string &inFile) {
	auto dump = load_dump(inFile);
	assert(!dump.empty());

	int skip_n_lines = 0;
	double total_weight = 0.;

	for (auto &line : dump) {
		if (skip_n_lines > 0){
			//it was set, so decrement and skip.
			skip_n_lines--;
			if (debug) fprintf(stderr,"RTPLan: skipping line: %s\n", line.second);
			continue;
		}
		if (line.first == "isocentername"){
			metaData.isocentername = line.second;
			continue;
		}
		if (line.first == "machinenameandversion") {
			if (startswith(line.second, "MLC160")){
				metaData.accelerator = Accelerator(AcceleratorType::Agility);
				continue;
			}
			else if (startswith(line.second, "M160")){
				metaData.accelerator = Accelerator(AcceleratorType::Agility); //FFF?
				continue;
			}
			else if (startswith(line.second, "MLC80")){
				metaData.accelerator = Accelerator(AcceleratorType::MLCi80);
				continue;
			}
			else if (startswith(line.second, "M80")){
				metaData.accelerator = Accelerator(AcceleratorType::MLCi80);
				continue;
			}
			else {
				metaData.accelerator = Accelerator(AcceleratorType::UNKNOWN);
				continue;
			}
		}
		if (line.first == "machineenergyname") {
			if (startswith(line.second, "6")){
				metaData.accelerator.energy = Energy::MV6;
			}
			else if (startswith(line.second, "10")){
				metaData.accelerator.energy = Energy::MV10;
			}
			else if (startswith(line.second, "7")){ //cant happen I guess
				metaData.accelerator.energy = Energy::MV7;
			}
			else {
				metaData.accelerator.energy = Energy::UNKNOWN;
				metaData.accelerator.filter = Filter::UNKNOWN;
			}
			if (endswith(line.second, "FFF")){
				metaData.accelerator.filter = Filter::NoFF;
				//TODO: imrtfilter="Compensator" == FFF?
			}
			continue;
		}
		if (line.first == "numberofcontrolpoints") {
			controlPoints.resize(stoi(line.second));
			for (int i = 0; i < stoi(line.second); i++){
				controlPoints[i].collimator.mlc.leftLeaves.resize(metaData.accelerator.leafs_per_bank); // leaf_per_bank should be set before
				controlPoints[i].collimator.mlc.rightLeaves.resize(metaData.accelerator.leafs_per_bank);
			}
			continue;
		}
		if (startswith(line.first, "gantry[")) {
			controlPoints[get_index(line.first)].beamInfo.gantryAngle = { stof(line.second), stof(line.second) };
			continue;
		}
		if (startswith(line.first, "couch[")) {
			controlPoints[get_index(line.first)].beamInfo.couchAngle = { 360-stof(line.second), 360-stof(line.second) };
			// we doen 360 - waarde, want Pinnacle...
			continue;
		}
		if (startswith(line.first, "collimator[")) {
			controlPoints[get_index(line.first)].beamInfo.collimatorAngle = { stof(line.second), stof(line.second) };
			continue;
		}
		if (startswith(line.first, "setbeamtype")) {
			if (startswith(line.second, "Dynamic Arc")){
				metaData.beamtype = BeamType::VMAT;
			}
			else if (startswith(line.second, "Step & Shoot MLC")){
				metaData.beamtype = BeamType::IMRT;
			}
			else if (startswith(line.second, "Static")){
				metaData.beamtype = BeamType::IMRT;
			}
			else {
				metaData.beamtype = BeamType::UNKNOWN;
			}
			continue;
		}
		if (startswith(line.first, "numberofpoints[")) {
			//use as check.
			if (stoi(line.second) == metaData.accelerator.leafs_per_bank){
				continue; //all ok
			}
			else {
				fprintf(stderr,"RTPLan: Invalid CPI detected: incorrect nr leafs: %s. Skipping CPI...\n", line.second);
				skip_n_lines = stoi(line.second) * 2; //skip this nr*2 (pairs) of lines
				continue;
			}
		}
		if (startswith(line.first, "points_element[")) {
			//loop backwards!
			int cpi = get_index(line.first);
			int leafindex = get_index(line.first, 2);
			float pos = stof(line.second);
			//fprintf(stderr,"points element %i %i %f\n", cpi, leafindex, pos);
			if (leafindex % 2 == 0){ //even, dus leftbank
				//std::cout << "left[" << cpi << "][" << leafindex / 2 << "] " << -pos;
				controlPoints[cpi].collimator.mlc.leftLeaves[metaData.accelerator.leafs_per_bank - 1 - (leafindex / 2)] = { -pos, -pos }; //minus want vanaf het midden!
			}
			else{
				//std::cout << "right[" << cpi << "][" << (leafindex - 1) / 2 << "] " << pos;
				controlPoints[cpi].collimator.mlc.rightLeaves[metaData.accelerator.leafs_per_bank - 1 - ((leafindex - 1) / 2)] = { pos, pos };
			}
			continue;
		}
		if (startswith(line.first, "leftjawposition[")) {
			controlPoints[get_index(line.first)].collimator.parallelJaw.j1 = { -stof(line.second), -stof(line.second) }; //J1 always most negative coord according to doc.
			controlPoints[get_index(line.first)].beamInfo.fieldMin.first = -stof(line.second) - metaData.fieldMargin; // pinnacle gives distance to center of colli, so add minus
			continue;
		}
		if (startswith(line.first, "rightjawposition[")) {
			controlPoints[get_index(line.first)].collimator.parallelJaw.j2 = { stof(line.second), stof(line.second) };
			controlPoints[get_index(line.first)].beamInfo.fieldMax.first = stof(line.second) + metaData.fieldMargin;
			continue;
		}
		if (startswith(line.first, "topjawposition[")) { //TRF reader: y axis swapped ? increases downward
			controlPoints[get_index(line.first)].collimator.perpendicularJaw.j2 = { stof(line.second), stof(line.second) }; //J1 always most negative coord according to doc.
			controlPoints[get_index(line.first)].beamInfo.fieldMax.second = stof(line.second) + metaData.fieldMargin;
			continue;
		}
		if (startswith(line.first, "bottomjawposition[")) {
			controlPoints[get_index(line.first)].collimator.perpendicularJaw.j1 = { -stof(line.second), -stof(line.second) };
			controlPoints[get_index(line.first)].beamInfo.fieldMin.second = -stof(line.second) - metaData.fieldMargin;
			continue;
		}
		if (startswith(line.first, "weight[")) { //weight sums to one per beam.
			controlPoints[get_index(line.first)].beamInfo.relativeWeight = stof(line.second); //relweight is relative to computation, so total should be 1, and per beam it is.
			total_weight += stof(line.second);
			continue;
		}
		if (startswith(line.first, "weight") && !startswith(line.first, "weight[")) { // beam weight, one per beam. summed over beams should be one.
			metaData.weight = stof(line.second);
			continue;
		}
	}

	if (debug) fprintf(stderr,"RTPLan: One beam loaded with a beam weight of %f, and cps weights sum to %f.\n", metaData.weight, total_weight);
}


void Parser::setCPIs() {
	if (metaData.accelerator.type == AcceleratorType::EMPTY || metaData.accelerator.type == AcceleratorType::UNKNOWN) throw std::pair<int, string>(27, "Undefined accelerator encountered.");
	if (metaData.accelerator.energy == Energy::UNKNOWN) throw std::pair<int, string>(21, "Unknown energy encountered.");
	if (metaData.accelerator.filter == Filter::UNKNOWN) throw std::pair<int, string>(23, "Unknown filter encountered.");
	if (metaData.beamtype == BeamType::UNKNOWN) throw std::pair<int, string>(22, "Unknown beamtype encountered.");

	metaData.prescriptiondose *= (metaData.dose_per_fraction ? 1 : metaData.nr_fractions);

	double total_weight = 0.f;
	// zet orientaties, isoc, correct weights
	for (auto &cp : controlPoints){
		cp.beamInfo.isoCenter.x = metaData.isoc.x;
		cp.beamInfo.isoCenter.y = metaData.isoc.y;
		cp.beamInfo.isoCenter.z = metaData.isoc.z;
		cp.collimator.mlc.orientation = ModifierOrientation::IECY;
		cp.collimator.perpendicularJaw.orientation = ModifierOrientation::IECX;

		if (metaData.accelerator.type == AcceleratorType::MLCi80){
			cp.collimator.parallelJaw.orientation = ModifierOrientation::IECY;
		}
		else {
			cp.collimator.parallelJaw.orientation = ModifierOrientation::NOT_PRESENT;
		}
		total_weight += cp.beamInfo.relativeWeight;
		//(beam)weight given in perc, to factor (/100), and then from Gy to cGy (*100)
		cp.beamInfo.relativeWeight *= metaData.weight * metaData.mu_per_fraction * (metaData.dose_per_fraction ? 1 : metaData.nr_fractions);
		//								^-beam weight		^- mu/frac						^- times nb frac or not.
	}
	assert((total_weight >= 0.9999) && (total_weight <= 1.0001)); //cp weigths within beam must sum to 1. only machine error is tolerated.

	//VMAT to dynamic arcs
	if (metaData.beamtype == BeamType::VMAT && !metaData.pinnacle_vmat_interpolation){
		for (int i = 0; i < (num_cps() - 1); i++){ //skip last CPI
			// van N CPIs naar N-1 segmenten.
			// kijk ACHTERUIT: skippen dus N=0. Want, laatste CP heeft weight nul
			// Let op: segment fieldMin,fieldMax moeten omhullende van naastgelegen CPIs worden

			controlPoints[i].beamInfo.collimatorAngle.second = controlPoints[i + 1].beamInfo.collimatorAngle.first; //first is already correct
			controlPoints[i].beamInfo.couchAngle.second = controlPoints[i + 1].beamInfo.couchAngle.first;
			controlPoints[i].beamInfo.gantryAngle.second = controlPoints[i + 1].beamInfo.gantryAngle.first;

			//fieldMin/Max is static!!!, take most outward of the two CPIs.
			controlPoints[i].beamInfo.fieldMin.first = std::min(controlPoints[i].beamInfo.fieldMin.first, controlPoints[i + 1].beamInfo.fieldMin.first);
			controlPoints[i].beamInfo.fieldMin.second = std::min(controlPoints[i].beamInfo.fieldMin.second, controlPoints[i + 1].beamInfo.fieldMin.second);
			controlPoints[i].beamInfo.fieldMax.first = std::max(controlPoints[i].beamInfo.fieldMax.first, controlPoints[i + 1].beamInfo.fieldMax.first);
			controlPoints[i].beamInfo.fieldMax.second = std::max(controlPoints[i].beamInfo.fieldMax.second, controlPoints[i + 1].beamInfo.fieldMax.second);

			//j1 == negative, so minimum
			controlPoints[i].collimator.perpendicularJaw.j1.second = controlPoints[i + 1].collimator.perpendicularJaw.j1.first;
			controlPoints[i].collimator.perpendicularJaw.j2.second = controlPoints[i + 1].collimator.perpendicularJaw.j2.first;
			controlPoints[i].collimator.parallelJaw.j1.second = controlPoints[i + 1].collimator.parallelJaw.j1.first;
			controlPoints[i].collimator.parallelJaw.j2.second = controlPoints[i + 1].collimator.parallelJaw.j2.first;

			for (int j = 0; j < metaData.accelerator.leafs_per_bank; j++){
				controlPoints[i].collimator.mlc.leftLeaves[j].second = controlPoints[i + 1].collimator.mlc.leftLeaves[j].first;
				controlPoints[i].collimator.mlc.rightLeaves[j].second = controlPoints[i].collimator.mlc.rightLeaves[j].first;
			}
		}
		controlPoints.pop_back(); //we moved every CP forward, so last one is now superfluous
	}
	//VMAT pinnacle interpol. no dynamic thingies, but only reweighting
	else if (metaData.beamtype == BeamType::VMAT && metaData.pinnacle_vmat_interpolation){

		for (int i = 0; i < (num_cps() - 1); i++){ // loop over CPIs except last
			// we take average weight of i and i + 1, because of the pinnacle shifted weights wrt to dicom and tabel 5 in NLCOMvStralingsdosimetrie Rapport 26
			// i=1 has weight zero to we get half weight of next CP anyway.
			controlPoints[i].beamInfo.relativeWeight = controlPoints[i].beamInfo.relativeWeight / 2. + controlPoints[i + 1].beamInfo.relativeWeight / 2.;
		}

		controlPoints.pop_back(); //weights have shifted forward, laatste mag weg
	}

}


void dicomParser::setCPIs() {
	//first, convert CUMULATIVE mu to mu.
	// kijk vooruit: skippen dus N=N. Want, eerste CP heeft weight 0
	// skip last, because lookahead on next line
	for (int i = 0; i < num_cps() - 1; i++){
		controlPoints[i].beamInfo.relativeWeight = controlPoints[i + 1].beamInfo.relativeWeight - controlPoints[i].beamInfo.relativeWeight; //i+1 - i
	}
	controlPoints[num_cps() - 1].beamInfo.relativeWeight = 0; //last one can be zeroed out.
	Parser::setCPIs(); //then regular procedure
}


//RTBeam::RTBeam(const string &rt_files, float _fieldMargin, bool _debug, bool _pinnacleVMATmode) {
RTBeam::RTBeam(DosiaSettings &_sett) : sett(_sett){
	string &rt_files = sett.rt_files;
	io::isfile(rt_files + "/dbtype.dump", 20);
	auto dbtype = load_dump(sett.rt_files + "/dbtype.dump");
	for (auto &line : dbtype) {
		if (line.first == "pinnacle"){
			auto parser = pinnacleParser(metaData, (sett.verbose > 1)?true:false);
			parser.parseTrial(rt_files + "/trialname.dump");	//mu, nr_fracties, HU-intercept. NIET beschikbaar bij dicom
			parser.parseBeam(rt_files + "/beam.dump");
			parser.parsePlan(rt_files + "/plan.dump");	//isoc, requires info from beam!
			parser.parseScan(rt_files + "/scan.dump");	//not needed
			parser.parseDose(rt_files + "/dose.dump");

			parser.setCPIs();
			metaData = parser.metaData;
			controlPoints = parser.controlPoints;
		}
		else if (line.first == "dicom"){
			throw std::pair<int, string>(28, "Dicom not fully implemented.");
			auto parser = dicomParser(metaData, (sett.verbose > 1) ? true : false);
			parser.parseTrial(rt_files + "/trialname.dump");	//NIET beschikbaar bij dicom
			parser.parsePlan(rt_files + "/plan.dump");	//no required info?
			parser.parseBeam(rt_files + "/beam.dump");
			parser.parseScan(rt_files + "/scan.dump");	//dicom slope/intercept, ImagePositionPatient
			parser.parseDose(rt_files + "/dose.dump"); //not needed

			parser.setCPIs();
			metaData = parser.metaData;
			controlPoints = parser.controlPoints;
		}
	}

	if (sett.verbose > 2) {
		printInfo();
		printFirstLeaf();
	}
};


void RTBeam::printInfo(){
	fprintf(stderr,"RTPLan: Number of fractions is %i, %.2f MU per fraction and a prescription dose of %.2f.\n", metaData.nr_fractions, metaData.mu_per_fraction, metaData.prescriptiondose);
	fprintf(stderr,"RTPLan: This beam of weight %.2f has %i controlpoints. Per CPI:\n", metaData.weight, num_cps());
	for (auto &cp : controlPoints){
		fprintf(stderr,"fieldMin/Max                 \t %.2f \t %.2f \t %.2f \t %.2f \n", cp.beamInfo.fieldMin.first, cp.beamInfo.fieldMin.second, cp.beamInfo.fieldMax.first, cp.beamInfo.fieldMax.second);
		fprintf(stderr,"parallelJaw.j1/j2/orient     \t %.2f \t %.2f \t %i \n", cp.collimator.parallelJaw.j1.first, cp.collimator.parallelJaw.j2.first, cp.collimator.parallelJaw.orientation);
		fprintf(stderr,"perpendicularJaw.j1/j2/orient\t %.2f \t %.2f \t %i \n", cp.collimator.perpendicularJaw.j1.first, cp.collimator.perpendicularJaw.j2.first, cp.collimator.perpendicularJaw.orientation);
		fprintf(stderr,"gantryangle %.2f\n", cp.beamInfo.gantryAngle.first);
		fprintf(stderr,"relativeWeight (in MU) %.2f\n", cp.beamInfo.relativeWeight);
	}
}


void RTBeam::printFirstLeaf(){
	fprintf(stderr,"RTPLan: This beam has %i controlpoints. First CPI:\n", num_cps());
	for (int i = 0; i < controlPoints[0].collimator.mlc.leftLeaves.size(); i++){
		fprintf(stderr,"%i th leafbank (l, r): %.2f, %.2f\n", i, controlPoints[0].collimator.mlc.leftLeaves[i].first, controlPoints[0].collimator.mlc.rightLeaves[i].first);
	}
}
