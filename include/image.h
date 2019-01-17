#pragma once

#include <iostream>
#include "pystring.h" //pystring
#include "tools.h" //vect,std::vector
#include <cstring> //std::memcpy
#include "gpumcd/Phantom.h"
using namespace vect;
using namespace pystring;

class Image{
public:
	vector<int> dim_size;
	vector<float> voxel_sizes;
	vector<float> min_ext;
	vector<float> max_ext;
	vector<float> imdata; // watch out! I convert all to floats! my world is simple!

	Image() = default;
	~Image() = default;
	Image(const std::string &);

	void write(const std::string &);
	Image copy_with_new_voxels(const vector<float> &);

	int ndim() const { return dim_size.size(); };
	int nvox() const { return mul(dim_size); };

private:
	void read_xdr(const std::string &);
	void read_mhd(const std::string &);

	void write_xdr(const std::string &);
	void write_mhd(const std::string &);
};


Image::Image(const std::string &fname){
	if (pystring::endswith(fname, ".xdr")){
		read_xdr(fname);
	}
	else if (pystring::endswith(fname, ".mhd")){
		read_mhd(fname);
	}
	else {
		printf("Unkown file extension encountered, aborting write out...");
	}
}


void Image::write(const std::string &fname){
	if (pystring::endswith(fname, ".xdr")){
		write_xdr(fname);
	}
	else if (pystring::endswith(fname, ".mhd")){
		write_mhd(fname);
	}
	else {
		printf("Unkown file extension encountered, aborting write out...");
	}
}


void Image::read_xdr(const std::string &xdrfile) {
	//Phantom phantom;

	auto xdr = fromfile<char>(xdrfile);

	std::string header;

	//calc offset of magic bytes
	char lasti = ' ';
	long imdata_offset = 0;
	for (char i : xdr){
		imdata_offset++;
		if (i == 0x0c && lasti == 0x0c){
			break;
		}
		lasti = i;
		header += i;
	}
	header.pop_back(); //one magic byte was added.

	int type = -1; //2 = >i2, 4 = >f4

	for (const auto &line : parse::parse_dump(header)){
		if (startswith(line.first, "ndim")) {
			int _ndim = stoi(line.second);
			assert(_ndim > 0 && _ndim < 4);
			dim_size.resize(_ndim);
			voxel_sizes.resize(_ndim);
			min_ext.resize(_ndim);
			max_ext.resize(_ndim);
			continue;
		}
		if (startswith(line.first, "field")) {
			assert(startswith(line.second, "uniform"));
			continue;
		}
		if (startswith(line.first, "data")) {
			if (startswith(line.second, "xdr_short")){
				type = 2;
			}
			else if (startswith(line.second, "xdr_real") || startswith(line.second, "xdr_float")){
				type = 4;
			}
			else{
				assert(type != -1); //blow up
			}
			continue;
		}
		if (startswith(line.first, "dim1")) {
			dim_size[0] = stoi(line.second);
			continue;
		}
		if (startswith(line.first, "dim2")) {
			dim_size[1] = stoi(line.second);
			continue;
		}
		if (startswith(line.first, "dim3")) {
			dim_size[2] = stoi(line.second);
			continue;
		}
		/*if (startswith(line.first, "min_ext")) {
			min_ext.push_back(stof(split(line.second)[0]));
			min_ext.push_back(stof(split(line.second)[1]));
			min_ext.push_back(stof(split(line.second)[2]));
			continue;
			}
			if (startswith(line.first, "max_ext")) {
			max_ext.push_back(stof(split(line.second)[0]));
			max_ext.push_back(stof(split(line.second)[1]));
			max_ext.push_back(stof(split(line.second)[2]));
			continue;
			}*/
	}

	//header loaded.

	long imdata_bytes = nvox()*type;
	long ext_offset = imdata_offset + imdata_bytes;
	long ext_bytes = ndim() * 2 * sizeof(float);

	//check that the size of the .xdr corresponds to the header+voxels*voxeltype+exts:
	assert(xdr.size() == imdata_offset + imdata_bytes + ext_bytes);

	vector<char> _imdata(xdr.begin() + imdata_offset, xdr.begin() + imdata_offset + imdata_bytes);

	if (type == 2){
		types::swap_endianness<short>(_imdata);
		vector<short> shorts = types::reinterpret<short>(_imdata);
		imdata.insert(imdata.begin(), shorts.begin(), shorts.end()); //this upcasts shorts
	}
	else if (type == 4){
		types::swap_endianness<float>(_imdata);
		imdata = types::reinterpret<float>(_imdata);
	}

	//now the extents in the final ndim*2*sizeof(float) bytes
	//write extents, looped pairwise over axis
	//xmin, xmax, ymin, ymax, zmin, zmax

	vector<char> _exts(xdr.begin() + ext_offset, xdr.begin() + ext_offset + ext_bytes);
	types::swap_endianness<float>(_exts);
	vector<float> exts = types::reinterpret<float>(_exts);

	for (size_t i = 0; i < ndim(); i++) {
		min_ext[i] = exts[2 * i];
		max_ext[i] = exts[2 * i + 1];
		//calc binsize
		voxel_sizes[i] = (max_ext[i] - min_ext[i]) / (dim_size[i] - 1);
	}
}


void Image::write_mhd(const std::string &fn){
	FILE* ffile = fopen(fn.c_str(), "wb");
	if (ffile == nullptr){
		throw std::pair<int, std::string>(70, "Problem writing file '" + fn + "'.");
	}
	// TODO:header
	fprintf(ffile, "ObjectType = Image\n");
	fprintf(ffile, "NDims=%i\n", ndim());
	fprintf(ffile, "BinaryData = True\n");
	fprintf(ffile, "BinaryDataByteOrderMSB = False\n");
	fprintf(ffile, "CompressedData = False\n");
	fprintf(ffile, "Offset = ");
	for (int i = 0; i < ndim(); i++) {
		fprintf(ffile, "%f ", min_ext[i]*10);
	}
	fprintf(ffile, "\n");
	fprintf(ffile, "ElementSpacing = ");
	for (int i = 0; i < ndim(); i++) {
		fprintf(ffile, "%f ", voxel_sizes[i]*10);
	}
	fprintf(ffile, "\n");
	fprintf(ffile, "DimSize = ");
	for (int i = 0; i < ndim(); i++) {
		fprintf(ffile, "%i ", dim_size[i]);
	}
	fprintf(ffile, "\n");
	fprintf(ffile, "ElementType = MET_FLOAT\n");
	std::string rawfile;
	std::string ext;
	os::path::splitext(rawfile,ext,fn);
	rawfile += ".raw";
	fprintf(ffile, "ElementDataFile = %s\n",rawfile.c_str());
	fclose(ffile);

	//write rawfile
	tofile<float>(imdata,rawfile);
}


void Image::write_xdr(const std::string &fn){
    FILE* ffile = fopen(fn.c_str(), "wb");
    if (ffile == nullptr){
        throw std::pair<int, std::string>(70, "Problem writing file '" + fn + "'.");
    }
    fprintf(ffile, "# AVS WRITER BY BRENT\n");
    fprintf(ffile, "ndim=%i\n",ndim());
	for (int i = 0; i < ndim(); i++) {
        fprintf(ffile, "dim%i=%i\n",i+1,dim_size[i]);
    }
    fprintf(ffile, "nspace=%i\n",ndim()); //dont know if used
    fprintf(ffile, "veclen=1\n"); //dont know if used
    fprintf(ffile, "data=xdr_real\n"); //dont know if used
    fprintf(ffile, "field=uniform\n"); //dont know if used
    fprintf(ffile, "%c%c",0x0c,0x0c); //magic bytes

    //convert imdata to bytes
	vector<char> imdata_swapped = types::reinterpret<char>(imdata);
	types::swap_endianness<float>(imdata_swapped);
    fwrite(imdata_swapped.data(), sizeof(char), imdata_swapped.size(), ffile);

	//extents
	//xmin, xmax, ymin, ymax, zmin, zmax
	vector<float> exts(ndim() * 2);
	for (size_t i = 0; i < ndim(); i++) {
		exts[2 * i]=min_ext[i];
		exts[2 * i + 1]=max_ext[i];
	}
	vector<char> exts_swapped = types::reinterpret<char>(exts);
	types::swap_endianness<float>(exts_swapped);
    fwrite(exts_swapped.data(), sizeof(char), exts_swapped.size(), ffile);
    fclose(ffile);
}


void Image::read_mhd(const std::string &header) {
	std::string rawfile;
	int type = -1; //2 = >i2, 4 = >f4
	for (const auto &line : parse::load_dump(header)) {
		if (startswith(line.first, "NDims")) {
			int _ndim = stoi(line.second);
			assert(_ndim > 0 && _ndim < 4);
			dim_size.resize(_ndim);
			voxel_sizes.resize(_ndim);
			min_ext.resize(_ndim);
			max_ext.resize(_ndim);
			continue;
		}
		if (startswith(line.first, "BinaryData") && !startswith(line.first, "BinaryDataByteOrderMSB")) {
			assert(startswith(line.second, "True"));
			continue;
		}
		if (startswith(line.first, "BinaryDataByteOrderMSB")) {
			assert(startswith(line.second, "False"));
			continue;
		}
		if (startswith(line.first, "CompressedData")) {
			assert(startswith(line.second, "False"));
			continue;
		}
		if (startswith(line.first, "ElementSpacing")) {
			voxel_sizes = types::split<float>(line.second);
			for (auto &i : voxel_sizes){
				i /= 10;
			}
			continue;
		}
		if (startswith(line.first, "DimSize")) {
			dim_size = types::split<int>(line.second);
			continue;
		}
		if (startswith(line.first, "Offset")) {
			min_ext = types::split<float>(line.second);
			for (auto &i : min_ext){
				i /= 10;
			}
			continue;
		}
		if (startswith(line.first, "ElementDataFile")) {
			rawfile = line.second;
			continue;
		}
		if (startswith(line.first, "ElementType")) {
			if (startswith(line.second, "MET_SHORT")){
				type = 2;
			}
			else if (startswith(line.second, "MET_FLOAT")){
				type = 4;
			}
			else{
				assert(type != -1); //blow up
			}
			continue;
		}
	}

	if (type == 2){
		vector<short> shorts = fromfile<short>(rawfile);
		imdata.insert(imdata.begin(), shorts.begin(), shorts.end());
	}
	else if (type == 4){
		imdata = fromfile<float>(rawfile);
	}
	assert(imdata.size() == nvox());

	for (size_t i = 0; i < ndim(); i++) {
		max_ext[i] = min_ext[i] + voxel_sizes[i] * (dim_size[i] -1);
	}
}