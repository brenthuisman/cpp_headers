#pragma once

//#include <set>
#include <assert.h>

#include "CalculationInformation.h"

#include <string>
#include "pystring.h"
using namespace pystring;
using std::string;
#include "tools.h"
using namespace parse;
using namespace vect;
using std::vector;
#include "settings.h"

#include "image.h"

class FluenceImage;

class FluenceBeam;


class FluenceImage : private Image {
	FluenceImage(const std::string &);// : Image(const string &);
	std::vector<BeamInformation> beaminfos;
};

FluenceImage::FluenceImage(const string &fn) : Image(fn) {
	assert(ndim() == 2);
	beaminfos = std::vector<BeamInformation>(nvox());

	//x-coord panel ligt tussen min_ext[0]+voxel_sizes[0]*x waar x<max_ext[0]
	/*float relativeWeight;
	Float3 isoCenter;
	std::pair<float, float> gantryAngle;
	std::pair<float, float> couchAngle;
	std::pair<float, float> collimatorAngle;
	std::pair<float, float> fieldMin;
	std::pair<float, float> fieldMax;*/

	for (int y = 0; y < dim_size[1]; y++){
		for (int x = 0; x < dim_size[0]; x++){
			// todo convert x,y afmetingen in isoc
			// plus and minus halfpixel
			beaminfos[x + y*x].relativeWeight = imdata[x + y*x];
			beaminfos[x + y*x].fieldMin = { min_ext[0] + (x - 0.5)*voxel_sizes[0], min_ext[1] + (y - 0.5)*voxel_sizes[1] };
			beaminfos[x + y*x].fieldMax = { min_ext[0] + (x + 0.5)*voxel_sizes[0], min_ext[1] + (y + 0.5)*voxel_sizes[1] };
		}
	}

};