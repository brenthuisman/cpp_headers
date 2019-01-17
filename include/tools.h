#pragma once

#include <iostream>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element, std::sort, std::reverse, std::transform
#include <functional> //std::plus, std::minus
#include <numeric>	// std::iota
#include <random>
#include <fstream>
#include <assert.h>
#include <sstream> //stringstream
#include <cstring> //std::memcpy
#include "pystring.h"
#pragma warning(disable : 4996) // disable fopen warning vs

namespace vect {
	using std::vector;

	template <typename T, typename U, typename V>
	double interpolate(const vector<T> &xData, const vector<U> &yData, const V &x, const bool extrapolate) {
		//Can call with any combo of input types, but result is always double.
		//extrapolate sets gradient outside xData to zero or continue from nearest bin.
		size_t size = xData.size();

		size_t i = 0;// find left end of interval for interpolation
		if (x >= xData[size - 2]) {// special case: beyond right end
			i = size - 2;
		}
		else {
			while (x > xData[i + 1]) i++;
		}
		double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1];
		if (!extrapolate) {	// if beyond ends of array and not extrapolating, set gradient to zero
			if (x < xL) yR = yL;
			if (x > xR) yL = yR;
		}
		double gradient = (yR - yL) / (xR - xL);

		return yL + gradient * (x - xL); // linear interpolation
	};

	template <typename T>
	vector<size_t> sort_indexes(const vector<T> &v) {
		// src: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes#12399290
		// as per std::sort convention, low to high. use std::reverse

		// initialize original index locations
		vector<size_t> idx(v.size());
		std::iota(idx.begin(), idx.end(), 0);

		// sort indexes based on comparing values in v
		std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

		return idx;
	}

	vector<double> uniform_random_vector(int length = 100) {
		std::random_device rng;
		std::uniform_real_distribution<double> distribution(0.0, 1.0);
		vector<double> result(length);
		for (auto& item : result) item = distribution(rng);
		return result;
	}

	template <typename T>
	T sum(vector<T> &v){
		T sum = 0;
		for (const auto &i : v) {
			sum += i;
		}
		return sum;
	}

	template <typename T>
	void reverse(vector<T> &v){ //REFERENCE, otherwise its copied!!! it reverses IN PLACE.
		std::reverse(v.begin(), v.end());
	}

	//element wise add two vectors
	template <typename T>
	vector<T> add(const vector<T> &v1, const vector<T> &v2){
		assert(v1.size() == v2.size());
		vector<T> v3(v1.size());
		std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(), std::plus<T>());
		return v3;
	}

	//element wise subtract two vectors
	template <typename T>
	vector<T> sub(const vector<T> &v1, const vector<T> &v2){
		assert(v1.size() == v2.size());
		vector<T> v3(v1.size());
		std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(), std::minus<T>());
		return v3;
	}

	//element wise multiply two vectors
	template <typename T>
	vector<T> mul(const vector<T> &v1, const vector<T> &v2){
		assert(v1.size() == v2.size());
		vector<T> v3(v1.size());
		std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(), std::multiplies<T>());
		return v3;
	}

	//multiply elements of vector
	template <typename T>
	T mul(const vector<T> &v1){
		T retval = 1;
		for (const auto &n : v1) {
			retval *= n;
		}
		return retval;
	}

	//element wise divide two vectors
	template <typename T>
	vector<T> div(const vector<T> &v1, const vector<T> &v2){
		assert(v1.size() == v2.size());
		vector<T> v3(v1.size());
		std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(), std::divides<T>());
		return v3;
	}

	//write vector to binary file
	template <typename T>
	void tofile(const vector<T> &v, const std::string &fn){
		FILE* ffile = fopen(fn.c_str(), "wb");
		if (ffile == nullptr){
			throw std::pair<int, std::string>(70, "Problem writing file '" + fn + "'.");
			//throw "Problem writing file '" + fn + "'.";
		}
		fwrite(reinterpret_cast<const char*>(&v[0]), sizeof(T), v.size(), ffile);
		fclose(ffile);
	}

	//template specialization for vector of std::strings
	template <>
	void tofile<std::string>(const vector<std::string> &v, const std::string &fn){
		FILE* ffile = fopen(fn.c_str(), "w");
		if (ffile == nullptr){
			throw std::pair<int, std::string>(71, "Problem writing file '" + fn + "'.");
			//throw "Problem writing file '" + fn + "'.";
		}
		for (const auto &line : v) {
			fprintf(ffile, "%s\n",line.c_str());
		}
		fclose(ffile);
	}

	//read vector from binary file, unknown size
	template <typename T>
	vector<T> fromfile(const std::string &fn){
		FILE* ffile = fopen(fn.c_str(), "rb");
		if (ffile == nullptr){
			throw std::pair<int, std::string>(72, "Problem reading file '" + fn + "'.");
		}

		//get filesize
		fseek(ffile, 0, SEEK_END);
		long fsize = ftell(ffile) / sizeof(T);
		rewind(ffile);

		vector<T> v(fsize);

		fread(&v[0], sizeof(T), v.size(), ffile);
		fclose(ffile);
		return v;
	}

	//fromfile specialization for textfile
	template <>
	vector<std::string> fromfile(const std::string &fn) {
		std::ifstream f(fn);// .c_str());
		if (!f.good()) throw std::pair<int, std::string>(73, "Problem reading file '" + fn + "'.");
		vector<std::string> ret;
		std::string str;
		while (std::getline(f, str)) {
			ret.push_back(str);
		}
		return ret;
	}

	//len()
	template <typename T>
	size_t len(const vector<T> &source) {
		return source.size();
	}

	//index()
	template <typename T>
	size_t index(const vector<T> &source, const T &item) {
		size_t index = std::distance(source.begin(), std::find(std::begin(source), std::end(source), item));
		if (index == len(source)) throw std::pair<int, std::string>(74, "Element not found");
		return index;
	}

	//printers voor (nested) vectors van any type
	void print(const vector<std::string> &v, std::ostream& dest = std::cout){
		for (const auto &i : v) {
			dest << i << std::endl;
		}
	}
	template <typename T>
	void print(const vector<vector<T>> &v, std::ostream& dest = std::cout){
		for (const auto &row : v) {
			for (const auto &col : row) {
				dest << col << "\t";
			}
			dest << std::endl;
		}
	}
	template <typename T>
	void print(const vector<T> &v, std::ostream& dest = std::cout){
		for (const auto &i : v) {
			dest << i << "\t";
		}
		dest << std::endl;
	}
	template <typename T>
	void print(const T &t, std::ostream& dest = std::cout){
		dest << t << std::endl;
	}
}

namespace io {
	bool isfile(const std::string& filename) {
		std::ifstream f(filename);// .c_str());
		return f.good();
	}

	bool isfile(const std::string& filename, int errcode) {
		if (!isfile(filename)) throw std::pair<int, std::string>(errcode, "Failed to load " + filename + ".");
		return true; //if no throw, then OK!
	}
}

namespace types {
	template <typename U,typename T>
	std::vector<U> reinterpret(const std::vector<T> &in){
		//change vector type but keep an exact copy of buffer
		std::vector<U> out;
		if (in.empty()) return out; //nothing to do
		size_t typesize = sizeof(T);
		size_t nbytes = in.size()*typesize;
		size_t out_nelems = nbytes / sizeof(U);
		//printf("%i %i %i %i\n",typesize,sizeof(U),nbytes,out_nelems);
		out.resize(out_nelems);
		std::memcpy(out.data(), in.data(), nbytes);
		return out;
	}

	template <typename T>
	void swap_endianness(std::vector<char> &in){ //in place
		if (in.empty()) return; //nothing to do
		size_t typesize = sizeof(T);
		if (typesize == 1) return; //no endianness with onebyte types

		std::vector<char> out(in);
		for (size_t i = 0; i < in.size(); i += typesize) {
			for (size_t j = 0; j < typesize; j++){
				in[i + j] = out[i + typesize - 1 - j];
			}
		}
		return;
	}

	// lexical_cast, default is fallback
	template<class T>
	T lexical_cast(const std::string &str)
	{
		static std::istringstream very_long_name_ss; /* reusing has severe (positive) impact on performance */
		T value;
		very_long_name_ss.str(str);
		very_long_name_ss >> value;
		very_long_name_ss.clear();
		return value;
	}

	// lexical_cast, trivial conversion
	template<> std::string lexical_cast(const std::string &str) { return str; }

	// lexical_cast, conversions that exist in stl
	template<> float lexical_cast(const std::string &str) { return std::stof(str); }
	template<> double lexical_cast(const std::string &str) { return std::stod(str); }
	template<> long double lexical_cast(const std::string &str) { return std::stold(str); }
	template<> int lexical_cast(const std::string &str) { return std::stoi(str); }
	template<> long lexical_cast(const std::string &str) { return std::stol(str); }
	template<> long long lexical_cast(const std::string &str) { return std::stoll(str); }
	template<> unsigned long lexical_cast(const std::string &str) { return std::stoul(str); }
	template<> unsigned long long lexical_cast(const std::string &str) { return std::stoull(str); }

	// cast, conversions that need to be truncated
	template<> short lexical_cast(const std::string &str) { return static_cast<short>(lexical_cast<long>(str)); }
	template<> unsigned short lexical_cast(const std::string &str) { return static_cast<unsigned short>(lexical_cast<unsigned long>(str)); }
	template<> unsigned int lexical_cast(const std::string &str) { return static_cast<unsigned int>(lexical_cast<unsigned long>(str)); }

	// cast while splitting std::string
	template <typename T>
	std::vector<T> split(const std::string &source, const std::string &delim = ""){
		std::vector<std::string> s;
		if (delim.size() == 0){
			s = pystring::split(source);
		}
		else {
			s = pystring::split(source, delim);
		}
		std::vector<T> ret;
		for (const auto &i : s){
			ret.push_back(lexical_cast<T>(i));
		}
		return ret;
	}
}

namespace parse {

	int get_index(const std::string &str, const int &occurance){
		if (str.find("[", 1) != std::string::npos) { // there are brackets!
			int index_first_bracket = 1; //skip first char
			for (int i = 0; i < occurance; i++){
				index_first_bracket = str.find("[", index_first_bracket) + 1;
			}
			int index_second_bracket = str.find("]", index_first_bracket);
			return stoi(str.substr(index_first_bracket, index_second_bracket - index_first_bracket));
		}
		else {
			return -1;
		}

	}
	int get_index(const std::string &str){ return get_index(str, 1); }

	std::vector<std::pair<std::string, std::string>> parse_dump(const std::vector<std::string> &source) {
		std::vector<std::pair<std::string, std::string>> ret;
		for (const auto &item : source){
			if (pystring::split(item, "=").size() == 2) { //cant have pairs otherwise
				std::string key(pystring::strip(pystring::split(item, "=")[0]));
				if (pystring::startswith(key, "\"")) continue;
				std::string value(pystring::strip(pystring::split(item, "=")[1]));
				if (pystring::startswith(value, "\"") && value.length() > 2) {
					value = value.substr(1, value.size() - 2); //strips quotation marks
				}
				ret.push_back({ key, value });
			}
			else { // no two parts
				ret.push_back({ item, "" });
			}
		}
		return ret;
	}

	std::vector<std::pair<std::string, std::string>> parse_dump(const std::string &source) {
		return parse_dump(pystring::split(source, "\n"));
	}

	std::vector<std::pair<std::string, std::string>> load_dump(const std::string &dumpfile) {
		return parse_dump(vect::fromfile<std::string>(dumpfile));
	}

}
