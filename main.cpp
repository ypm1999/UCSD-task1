#include<vector>
#include<map>
#include<iostream>
#include<iterator>
#include<algorithm>
#include<chrono>
#include<fstream>
#include <sstream>
#include "threadPool.h"


#ifndef __linux__
#include "windows.h"
#else

#include "unistd.h"
#include "sys/sysinfo.h"
#endif

#ifndef __linux__
SYSTEM_INFO sysInfo;
GetSystemInfo( &sysInfo );
static const int CPU_NUM=sysInfo.dwNumberOfProcessors;
#else
static const int CPU_NUM=get_nprocs();
#endif
threadPool *threads = new threadPool(CPU_NUM);



typedef unsigned int SID;
typedef unsigned int Intensity;
typedef unsigned int MZ;
typedef unsigned int QID;
MZ const MZ_NONE=-1;

typedef std::pair<MZ, Intensity> Peak;
typedef std::vector<Peak> Spectrum;
typedef std::vector<Spectrum> RawData;

typedef std::pair<SID, Intensity> BucketPeak;

typedef std::vector<BucketPeak> Bucket;

typedef std::vector<Bucket> Index;
typedef std::tuple<> Answer;

static const MZ MAX_MZ = 20000;
static const int num_buckets = 1; //# of bins per mz used for index


class myAnswer{

public:
    MZ mz;
    QID qid;
    SID sid;
    Intensity intesity;
    myAnswer(){}

    myAnswer(MZ m, QID q, SID s, Intensity inten):
            mz(m), qid(q), sid(s), intesity(inten){}

    bool operator<(const myAnswer &b){
      if (mz == b.mz){
        return qid < b.qid;
      }
      else
        return mz < b.mz;

    }
};

RawData * load_raw_data(char *file, int &total_spectra, int &num_peaks) {
	RawData * spectra = new RawData();
	unsigned int file_id;
	unsigned int num_spectra;           //total number of spectra inside the file
	std::vector< std :: pair<unsigned int, unsigned int> > position; // starting position and ending position for each spectrum
	std::ifstream in(file, std::ios::in | std::ios::binary);
	in.read((char*)&file_id, sizeof( unsigned int ));
	in.read((char*)&num_spectra, sizeof( unsigned int ));
	position.resize(num_spectra);

	//total_spectra = 0 as input means to load all spectra from the DB
	if (total_spectra == 0){
		total_spectra = num_spectra;
	}

	if (num_spectra > 0) {
		position[0].first = num_spectra + 2; //starting offset in the file of the first spectrum
		in.read((char *) &position[0].second, sizeof(unsigned int)); //ending position of the first spectrum

		for (unsigned int spec_idx = 1; spec_idx < num_spectra; ++spec_idx) {
			position[spec_idx].first = position[spec_idx - 1].second;  //starting position
			in.read((char *) &position[spec_idx].second, sizeof(unsigned int)); //ending position
		}

		for (unsigned int spec_idx = 0; spec_idx < total_spectra; ++spec_idx) {
			in.seekg(position[spec_idx].first * 4 + 16); //setting the offset to read the spectrum
			unsigned int size = position[spec_idx].second - position[spec_idx].first - 4; //total number of peaks inside the spectrum

			//num_peaks = 0 as input means to load all peaks of the spectrum
			if (num_peaks == 0 || num_peaks < size){
				num_peaks = size;
			}

			Spectrum spectrum;
			//Populating peak info per spectrum
			for (int i = 0; i < num_peaks; ++i) {
				unsigned int peak;
				unsigned int mz;
				in.read((char *) &peak, sizeof(unsigned int));
				mz = peak >> 8;
				spectrum.push_back(Peak(mz,peak - (mz<<8)));
			}
			spectra -> push_back(spectrum);
		}
	}
	in.close();
	return spectra;
}

void dump_spectrum(Spectrum *s) {
	std::cerr << "[";
	for(auto & p: *s) {
		std::cerr << "[" << p.first << ", " << p.second << "],";
	}
	std::cerr << "]";
}

void dump_raw_data(RawData* r) {

	int c = 0;
	for(auto & s: *r) {
		std::cerr << "SID=" << c;
		dump_spectrum(&s);
		c++;
		std::cerr << "\n";
	}
}

void dump_index(Index *index) {
	for(MZ mz = 0; mz < MAX_MZ; mz++) {
		if(!(*index)[mz].empty()) {
			std::cerr << "MZ = " << mz << ": ";
			dump_spectrum(&(*index)[mz]);
			std::cerr << "\n";
		}
	}
}


void json_reconstruction(char * file, const std::vector<std::vector<myAnswer>> &reconstructed_spectra) {
	std::ofstream out(file, std::ios::out);
	
	out << "[\n";
	for(const auto & queries_results: reconstructed_spectra) {
		out << "[\n";
		MZ last_mz = MZ_NONE;
		for(const auto &answer:queries_results) {
		  if (last_mz != answer.mz){
		    if (last_mz != MZ_NONE){
          out << "\n";
          out << "\t\t]\n";
          out << "\t},\n";
		    }
        out << "\t{ " << "";
        out << "\"" << answer.mz << "\": " << "[\n";
		  }
		  else
		    out << ",\n";

      out << "\t\t\t[ " << answer.sid << ", " << answer.intesity << " ]";
      last_mz = answer.mz;
			}
    out << "\n";
    out << "\t\t]\n";
    out << "\t}\n";
		out << "]";
		if (&queries_results != &*reconstructed_spectra.rbegin()) {
			out << ",";
		}
		out << "\n";
	}
	out << "]\n";
}


Spectrum * load_query(char*file) {
	Spectrum *n = new Spectrum;
	std::ifstream in(file);
	std::string line;

	while (std::getline(in, line)) {
		std::vector<unsigned int> lineData;
		std::stringstream lineStream(line);
		unsigned int value;
		while (lineStream >> value) {
			lineData.push_back(value);
		}
		n->push_back(Peak(lineData[0], lineData[1]));
	}
	return n;
}

std::vector<Spectrum> * load_queries(char*file) {
	auto n = new std::vector<Spectrum>;

	std::ifstream in(file);
	std::string line;

	int spectra_count = 0;
	in >> spectra_count;
	for(int i = 0; i < spectra_count; i++) {
		Spectrum s;
		int peak_count;
		in >> peak_count;
		for (int p = 0; p < peak_count; p++) {
			MZ mz;
			Intensity intensity;
			in >> mz >> intensity;
			s.push_back(Peak(mz, intensity));
		}
		n->push_back(s);
	}
	return n;
}

Index * build_index(RawData * data) {
	Index *index = new Index;

	for(MZ mz = 0; mz < MAX_MZ; mz++) {
		index->push_back(Bucket());
	}

	unsigned int unit_frag;

	for(SID sid = 0; sid < data->size(); sid++) {
		for(auto & peak: (*data)[sid]) {
			unit_frag = peak.first/num_buckets;
			if (unit_frag < MAX_MZ) {
				(*index)[unit_frag].push_back(BucketPeak(sid, peak.second));
			}
		}
	}

	for(MZ mz = 0; mz < MAX_MZ; mz++) {
		std::sort((*index)[mz].begin(), (*index)[mz].end());
	}

	return index;
}


std::vector<myAnswer> reconstruct_candidate(Index * index, Spectrum query){
  std::vector<myAnswer> m;
  m.reserve(query.size());
  unsigned int cnt = 0;
  for(auto & query_peak: query) {
    unsigned int unit_q_mz = query_peak.first / num_buckets;
    for(auto & bucket_peak : (*index)[unit_q_mz]) {
      m.emplace_back(myAnswer(bucket_peak.first, cnt++, query_peak.first, bucket_peak.second));
    }
  }
  sort(m.begin(), m.end());
  return m;
}

std::vector<std::vector<myAnswer>> *reconstruct_candidates(Index * index, const std::vector<Spectrum> & queries) {
	auto reconstructed_spectra = new std::vector<std::vector<myAnswer>>;
  reconstructed_spectra->reserve(queries.size());
  std::queue< std::future<std::vector<myAnswer>> > results;
	for(auto & query: queries) {
	  if (results.size() > 4 * CPU_NUM){
      reconstructed_spectra->push_back(results.front().get());
      results.pop();
	  }
    results.push(threads->commit(reconstruct_candidate, index, query));
	}
	while (!results.empty()){
    reconstructed_spectra->push_back(results.front().get());
    results.pop();
  }

	return reconstructed_spectra;
}


int main(int argc, char * argv[]) {

	if (argc != 4) {
		std::cerr << "Usage: main <raw data file> <query file> <output json>\n";
		exit(1);
	}

	/* Declare number of spectra you want to load from the database
	 * 0 - load all spectra from the file
	 * n - load first N spectra from the file
	 * for initial test in example n = 3
	 */
	int total_spectra = 0;

	/* Declare number of peaks you want to load from each spectrum
	 * 0 - load all peaks from the spectrum
	 * m - load first N peaks from the spectrum
	 * for initial test in example m = 5
	 */
	int num_peaks = 0;

	bool demo = false;

	if (std::getenv("DEMO")) {
		total_spectra = 100000;
		num_peaks = 0;
		demo = true;
	}

	RawData * raw_data = load_raw_data(argv[1], total_spectra, num_peaks);
	//std::cerr << "raw_data=\n";
	//dump_raw_data(raw_data);

	std::vector<Spectrum> *queries = load_queries(argv[2]);
	if (demo) {
		std::cerr << "queries=\n";
		for(auto &s: *queries) {
			dump_spectrum(&s);
			std::cerr << "\n";
		}
		std::cerr << "\n";
	}

	// Here's where the interesting part starts

	auto index_build_start = std::chrono::high_resolution_clock::now();
	Index * index = build_index(raw_data);
	auto index_build_end = std::chrono::high_resolution_clock::now();

	delete raw_data;

//	if (demo) {
//		std::cerr << "Index: ";
//		dump_index(index);
//	}

	auto reconstruct_start = std::chrono::high_resolution_clock::now();
	auto reconstructed_spectra = reconstruct_candidates(index, *queries);
	auto reconstruct_end = std::chrono::high_resolution_clock::now();
	json_reconstruction(argv[3], *reconstructed_spectra);

	std::cerr << "Found " << reconstructed_spectra->size() << " candidates \n";
	std::cerr << "Building the index took " << (std::chrono::duration_cast<std::chrono::nanoseconds>(index_build_end - index_build_start).count()+0.0)/1e9 << " s\n";
	std::cerr << "Reconstruction took     " << (std::chrono::duration_cast<std::chrono::nanoseconds>(reconstruct_end - reconstruct_start).count()+0.0)/1e9 << " s\n";
}

