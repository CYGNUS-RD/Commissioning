// This reads a .hst midas history file and makes a plot. Useful to perform analysis on history data
// This is a template and can be further complicated reading multiple variables, multiple files
//Compile as g++ -std=c++17 Readhistory.cxx -o read_hst if you do not use ROOT for plotting (change code accordingly)
//Compile as g++ Readhistory.cxx -o read_hst `root-config --libs --cflags` if you use ROOT (it needs to be compiled with at least c++17)

#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <TApplication.h>
#include <TGraph.h>

constexpr uint32_t RT_DEF  = 0x46445348;
constexpr uint32_t RT_DATA = 0x41445348;

using namespace std;

struct HistRecord {
    uint32_t record_type;
    uint32_t event_id;
    uint32_t time;
    uint32_t def_offset;
    uint32_t data_size;
};

struct Tag {
    char name[32];
    uint32_t type;
    uint32_t n_data;
};

struct TagInfo {
    std::string name;
    uint32_t type;
    uint32_t ndata;
    uint32_t offset;
};

int tid_size[] = {
    0,1,1,1,2,2,4,4,4,4,8,1,0,0,0,0,0
};

int main(int argc, char** argv) {

    if(argc>4)
    {
        cerr<<"Error in using code!\nCorrect use: ./nameprog.exe <filepath_with_name> <event_id> <variable>\n(event_id is the equipment you want to read the history of, variable is the name of the feature to read)\n";
        return 1;
    }

    string filename = argv[1];
    int event_id_look = stoi(argv[2]);
    string variable = argv[3];

    std::ifstream f(filename.c_str(), std::ios::binary);
    if (!f) {
        cerr << "Cannot open file\n";
        return 1;
    }

    unordered_map<int,vector<TagInfo>> defs;

    vector<double> time;
    vector<double> y;

    HistRecord r;

    //Allows to open interactive window
    TApplication *app =  new TApplication("app",0,NULL);

    while (f.read(reinterpret_cast<char*>(&r), sizeof(r))) {

        if (r.record_type == RT_DEF) {

            char event_name[32];
            f.read(event_name, sizeof(event_name));

            int consumed = sizeof(event_name);
            int offset = 0;
            std::vector<TagInfo> tags;

            while (consumed < static_cast<int>(r.data_size)) 
            {

                Tag t;
                f.read(reinterpret_cast<char*>(&t), sizeof(t));

                tags.push_back({
                    std::string(t.name),
                    t.type,
                    t.n_data,
                    static_cast<uint32_t>(offset)
                });

                offset += t.n_data * tid_size[t.type];
                consumed += sizeof(Tag);
            }

            defs[r.event_id] = tags;
        }

        else if (r.record_type == RT_DATA) 
        {

            vector<char> buf(r.data_size);
            f.read(buf.data(), buf.size());

            if (r.event_id != event_id_look)
                continue;

            //cout << "\nTimestamp " << r.time << "\n";
            time.push_back(r.time);

            auto &tags = defs[event_id_look];

            for (auto &tag : tags) 
            {

                const char* p = buf.data() + tag.offset;

                //cout << tag.name << " = ";
                if(tag.name == variable)
                {
                    switch(tag.type){

                        case 6:
                            //std::cout << *reinterpret_cast<const uint32_t*>(p);
                            y.push_back(*reinterpret_cast<const uint32_t*>(p));
                            break;

                        case 7:
                            //std::cout << *reinterpret_cast<const int32_t*>(p);
                            y.push_back(*reinterpret_cast<const int32_t*>(p));
                            break;

                        case 9:
                            //std::cout << *reinterpret_cast<const float*>(p);
                            y.push_back(*reinterpret_cast<const float*>(p));
                            break;

                        case 10:
                            //std::cout << *reinterpret_cast<const double*>(p);
                            y.push_back(*reinterpret_cast<const double*>(p));
                            break;

                        default:
                            cout << "[type " << tag.type << "]";
                    }
                }

                //cout << "\n";
            }
        }

        else {
            f.seekg(r.data_size, std::ios::cur);
        }
    }
    auto g1 = make_shared<TGraph>(y.size(),time.data(),y.data());
    g1->SetTitle(Form("%s;UNIX t (s);Current (#muA)",variable.c_str()));
    g1->Draw("apl");

    app->Run(true);

}