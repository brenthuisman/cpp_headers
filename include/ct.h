#pragma once

#include <istream>
#include <assert.h>

#include "pystring.h"
using namespace pystring;
using std::string;
#include "tools.h"
using namespace vect;
#include "image.h"

#include "Phantom.h"
#include "settings.h"


/*
 * Sets up and returns a (GPUMCD) phantom.
 */


class CT;

class CT{
public:
	//members needed by gpumcd
	Phantom phantom;
	vector<string> materials;

	//ctors
	CT() = default;
	CT(const CT &) = default;
	~CT() = default;

	CT(DosiaSettings &, BeamMetaData &);

	//methods
	int num_vox(){ return image.nvox(); };
	Image generate_image(const vector<float> &);

private:
	//members
	BeamMetaData beamMetaData;
	DosiaSettings sett;
	Image image;

	//methods
	Phantom generate_phantom(const Image &, const string &, const string &);
	template <typename T>
	void set_hu2density(const string &, const vector<T> &, vector<float> &);
	//void set_hu2material(const string & = "hu2mat.ini"); //use with Schneider data from Gate
	void set_density2material(const string &, const vector<float> &, vector<float> &);
};


CT::CT(DosiaSettings &_sett, BeamMetaData &_beamMetaData) : beamMetaData(_beamMetaData), sett(_sett){

	string &rt_files = sett.rt_files;
	image = Image(os::path::join(rt_files, "ct.xdr"));
	if (sett.verbose > 1) cerr << "phantom file loaded: " << os::path::join(rt_files, "ct.xdr") << "\n";
	
	phantom = generate_phantom(image, os::path::join(sett.hounsfield_conversion_dir, "hu2dens.ini"), os::path::join(sett.hounsfield_conversion_dir, "dens2mat.ini"));

	if (sett.verbose > 1) {
		fprintf(stderr, "voxelSizes: %.2f,%.2f,%.2f\n", phantom.voxelSizes.x, phantom.voxelSizes.y, phantom.voxelSizes.z);
		fprintf(stderr, "numVoxels: %i,%i,%i\n", phantom.numVoxels.x, phantom.numVoxels.y, phantom.numVoxels.z);
		fprintf(stderr, "phantomCorner: %.2f,%.2f,%.2f\n", phantom.phantomCorner.x, phantom.phantomCorner.y, phantom.phantomCorner.z);
	}

	if (sett.dbgoutput){
		Image mediumIndex = generate_image(phantom.mediumIndexArray);
		Image massDensity = generate_image(phantom.massDensityArray);
		mediumIndex.write(os::path::join(rt_files, "mediumIndex.xdr"));
		massDensity.write(os::path::join(rt_files, "massDensityArray.xdr"));
	}
};


Phantom CT::generate_phantom(const Image &im, const string &hu2dens_fname, const string &dens2mat_fname){
	//think of this function as a constructor for Phantom structures
	assert(im.ndim() == 3);

	Phantom phantom;
	vector<float> ct_voxels; // gpumcd cares only about floats

	//setup coords,dimensions

	phantom.numVoxels.x = im.dim_size[0];
	phantom.numVoxels.y = im.dim_size[1];
	phantom.numVoxels.z = im.dim_size[2];

	phantom.phantomCorner.x = im.min_ext[0];
	phantom.phantomCorner.y = im.min_ext[1];
	phantom.phantomCorner.z = im.min_ext[2];

	phantom.voxelSizes.x = im.voxel_sizes[0];
	phantom.voxelSizes.y = im.voxel_sizes[1];
	phantom.voxelSizes.z = im.voxel_sizes[2];

	//now, move phantomCorner a half voxel

	phantom.phantomCorner.x -= phantom.voxelSizes.x / 2.;
	phantom.phantomCorner.y -= phantom.voxelSizes.y / 2.;
	phantom.phantomCorner.z -= phantom.voxelSizes.z / 2.;

	//convert image to HU units
	ct_voxels.resize(im.nvox());
	for (int i = 0; i < im.nvox(); i++){
		ct_voxels[i] = im.imdata[i] * beamMetaData.hu_slope + beamMetaData.hu_intercept;
	}
	//convert HU units to mass density and materials.
	set_hu2density(hu2dens_fname, ct_voxels, phantom.massDensityArray);
	set_density2material(dens2mat_fname, phantom.massDensityArray, phantom.mediumIndexArray);

	if (sett.in_aqua_vivo || sett.score_and_transport_in_water || sett.score_dose_to_water) {
		// set ref medium to water
		materials.push_back("Water"); //preserve existing materials, because GPUMCD still checks presence.
		sett.physicsSettings.referenceMedium = materials.size() - 1;
	}
	if (sett.in_aqua_vivo || sett.score_and_transport_in_water) {
		// set all voxels to water (we added it above at the end of 'materials')
		std::fill(phantom.mediumIndexArray.begin(), phantom.mediumIndexArray.end(), materials.size() - 1);
	}
	if (sett.in_aqua_vivo) {
		// set all voxels density 1. FIXME: currently sets all to 1, must set 1 only to patient mask
		
		//std::fill(phantom.massDensityArray.begin(), phantom.massDensityArray.end(), 1.f);
	}
	
	return phantom;
}


Image CT::generate_image(const vector<float> &new_voxels){
	assert(new_voxels.size() == image.nvox());
	Image ret(image);
	ret.imdata = new_voxels;
	//ret.imdata.insert(ret.imdata.begin(), new_voxels.begin(), new_voxels.end());
	return ret;
}


template <typename T>
void CT::set_hu2density(const string &fname, const vector<T> &ct_voxels, vector<float> &massDensityArray) {
	io::isfile(fname, 43);

	vector<float> density_hu;
	vector<float> density;

	std::ifstream is(fname);
	string str;
	while (getline(is, str))
	{
		split(str);
		density_hu.push_back( stoi(split(str)[0]) );
		density.push_back(stof(split(str)[1]));
	}

	massDensityArray.resize(num_vox());

	for (int i = 0; i < num_vox(); i++){
		massDensityArray[i] = interpolate(density_hu, density, ct_voxels[i], true);
		if (massDensityArray[i] < 0) massDensityArray[i] = 0.f;
	}
};


void CT::set_density2material(const string &fname, const vector<float> &massDensityArray, vector<float> &mediumIndexArray) {
	assert(massDensityArray.size()>0);//this ensures the densities are available.
	io::isfile(fname, 45);
	
	//use dens2mat with data from AvL clinic and Monaco defaults in appendix C of research manual.
	//NOTE: indices are continuous: a matindex of 1.4 will be a mix of 40% material 1 and 60% material 2 by weight.
	//gpumcd uses continuous material indices to mix materials. this vector helps.
	vector<int> continuous_material_index_axis;
	int i = 0;

	//load file into following. 'matireals' is a class member, because gpumcd later needs it
	vector<std::pair<int, int>> material_hu;
	vector<float> material_dens; //if used, indices must map to 'material' member

	std::ifstream is(fname);
	string str;
	while (getline(is, str))
	{
		split(str);
		material_dens.push_back(stof(split(str)[0]));
		materials.push_back(split(str)[1]);
		continuous_material_index_axis.push_back(i++);
	}

	mediumIndexArray.resize(num_vox());

	if (sett.continous_materials){
		for (int i = 0; i < num_vox(); i++){
			mediumIndexArray[i] = interpolate(material_dens, continuous_material_index_axis, massDensityArray[i], false);
			//if (phantom.mediumIndexArray[i] < 1) phantom.mediumIndexArray[i] = 0.f; //clip first materials (probably air)
		}
	}
	else { //integer materials
		for (int i = 0; i < num_vox(); i++){
			mediumIndexArray[i] = 0; //any densities lower than minimum found in file are assumed to be of material at first line.
			for (int j = 0; j < material_dens.size(); j++){
				if (material_dens[j]>massDensityArray[i]){
					break;
				}
				mediumIndexArray[i] = j;
			}
		}
	}
};



/*
void CT::set_hu2material(const string &fname){
	//not used, for compat with Gate
	io::isfile(fname, 44);

	vector<std::pair<int, int>> material_hu;
	vector<float> material_dens; //if used, indices must map to 'material' member


	std::ifstream is(fname);
	string str;
	while (getline(is, str))
	{
		split(str);
		material_hu.push_back({ stoi(split(str)[0]), stoi(split(str)[1]) });
		materials.push_back(split(str)[2]);
	}

	phantom.mediumIndexArray.resize(num_vox());

	for (int i = 0; i < num_vox(); i++){
		for (int j = 0; j < material_hu.size(); j++){
			if (material_hu[j].first < ct_voxels[i] && ct_voxels[i] < material_hu[j].second){
				phantom.mediumIndexArray[i] = j;
			}
		}
	}
};*/

