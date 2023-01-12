#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <sys/types.h>

#include "include/Area_maximization_minimization.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std::chrono;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;               // Point_2 object
typedef K::Segment_2 Segment_2;           // Segment_2 object
typedef CGAL::Polygon_2<K> Polygon_2;     // Polygon_2 object
typedef std::vector<Point_2> Points;      // vector with Point_2 objects
typedef std::vector<Segment_2> segments;  // vector with Segment_2 objects
typedef std::vector<Polygon_2> Polygon_v; // vector with Polygon_2 objects

typedef CGAL::Search_traits_2<K> T;
typedef CGAL::Fuzzy_iso_box<T> box;
typedef CGAL::Kd_tree<T> tree;

typedef std::vector<double> dist; // vector with distances from a point to an edge
typedef std::vector<int> areas;   // vector with polugon areas
typedef std::vector<int> findd;   // vector with position of visible edges from an interior point(position in polygon chain)
typedef std::vector<int> Areas;
typedef std::vector<Point_2>::iterator pveciterator;  // iterator gia vector apo points
typedef std::vector<Segment_2>::iterator segiterator; // iterator gia segments
typedef std::vector<double> distance;

int flag_algo = -1;    // algorithm we choose
int flag_min_max = -1; // minimization or maximixation
int option = -1;       // global local or subdivision
int L = -1;            // path length
int flagalgo = -1;
int flaginit = -1;
int flagedge = -1;
double threshold = -1; // threshold
Polygon_2 p;           // polygon
segments chain;        // chain
Points points;         // points

int main(int argc, char *argv[])
{
    srand((unsigned)time(NULL)); // for the random number i create for random edge selection
    int file_name_counter = 10;
    int cut_off;
    int how_many_points;

    std::string file_name_1 = "uniform-";
    std::string file_name_2 = "-1";
    std::string file_name_3 = "-2";
    std::string file_name_4 = ".instance";
    std::string file_name_5;
    std::string final_file_name;
    double ratio_array[5];
    std::ofstream project_out;
    project_out.open(argv[4]);
    project_out << "||\t<convex_hull random local_search>\t||\t<convex_hull min local_search>\t||\t<convex_hull max local_search>\t||";
    project_out << "\t<incremental random 1a local_search>\t||\t<incremental random 1b local_search>\t||\t<incremental random 2a local_search>\t||";
    project_out << "\t<incremental random 2b local_search>\t||\t<incremental min 1a local_search>\t||\t<incremental min 1b local_search>\t||";
    project_out << "\t<incremental min 2a local_search>\t||\t<incremental min 2b local_search>\t||"; // telos local search
    project_out << "\t<convex_hull random simulated_annealing local_search>\t||\t<convex_hull min simulated_annealing local_search>\t||" ;
    project_out << "\t<convex_hull max simulated_annealing local_search>\t||\t<convex_hull random simulated_annealing global_search>\t||";
    project_out << "\t<convex_hull min simulated_annealing global_search>\t||\t<cconvex_hull max simulated_annealing global_search>\t||" ;
    project_out << "\t<incremental random 1a simulated_annealing local_search>\t||\t<incremental random 1b simulated_annealing local_search>\t||" ;
    project_out << "\t<incremental random 2a simulated_annealing local_search>\t||\t<incremental random 2b simulated_annealing local_search>\t||" ;
    project_out << "\t<incremental min 1a simulated_annealing local_search>\t||\t<incremental min 1b simulated_annealing local_search>\t||" ;
    project_out << "\t<incremental min 2a simulated_annealing local_search>\t||\t<incremental min 2b simulated_annealing local_search>\t||" ;
    project_out << "\t<incremental max 1a simulated_annealing local_search>\t||\t<incremental max 1b simulated_annealing local_search>\t||" ;
    project_out << "\t<incremental max 2a simulated_annealing local_search>\t||\t<incremental max 2b simulated_annealing local_search>\t||" ;
    project_out << "\t<incremental random 1a simulated_annealing global_search>\t||\t<incremental random 1b simulated_annealing global_search>\t||" ;
    project_out << "\t<incremental random 2a simulated_annealing global_search>\t||\t<incremental random 2b simulated_annealing global_search>\t||" ;
    project_out << "\t<incremental min 1a simulated_annealing global_search>\t||\t<incremental min 1b simulated_annealing global_search>\t||" ;
    project_out << "\t<incremental min 2a simulated_annealing global_search>\t||\t<incremental min 2b simulated_annealing global_search>\t||" ;
    project_out << "\t<incremental max 1a simulated_annealing global_search>\t||\t<incremental max 1b simulated_annealing global_search>\t||" ;
    project_out << "\t<incremental max 2a simulated_annealing global_search>\t||\t<incremental max 2b simulated_annealing local_search>\t||" << std::endl;

    project_out << "Size||min_score\t||max_score\t||min_bound\t||max_bound\t||<<<<<<<<<<<<||" << std::endl;
    double *min;
    while (file_name_counter <= 100000)
    {
        project_out << std::endl <<file_name_counter << "\t||";
        if (file_name_counter < 100)
        {
            file_name_5 = "00000";
        }
        else if (file_name_counter < 1000)
        {
            file_name_5 = "0000";
        }
        else if (file_name_counter < 10000)
        {
            file_name_5 = "000";
        }
        else if (file_name_counter < 100000)
        {
            file_name_5 = "00";
        }
        else
        {
            file_name_5 = "0";
        }
        int r_number = 1 + (rand() % 2); // pick between instance 1 or 2
        if (r_number == 1)
        {
            std::string conv = argv[2];
            std::string s = std::to_string(file_name_counter);
            final_file_name = conv + file_name_1 + file_name_5 + s + file_name_2 + file_name_4;
            std::cout << final_file_name << std::endl;
        }
        else
        {
            std::string conv = argv[2];
            std::string s = std::to_string(file_name_counter);
            final_file_name = conv + file_name_1 + file_name_5 + s + file_name_3 + file_name_4;
            std::cout << final_file_name << std::endl;
        }
        std::ifstream in_f(final_file_name);
        cut_off = 500 * file_name_counter;
        std::string line;
        std::getline(in_f, line); // skip the first two lines from input(we dont need them)
        std::getline(in_f, line);
        Points result2;
        result2 = handleinput(in_f, result2); // function tha returns a vector with all the points given from the input file
        std::ofstream outfile;
        for (int j = 1; j <= 3; j++) // convex hull random,min,max local search min,max
        {
            outfile.open("output_polygon.txt", std::ofstream::out | std::ofstream::trunc); // output file
            flagedge = j;                                                                  // random ,min ,max
            int min_m = 1;                                                                 // min
            convex_hull_fun(result2, outfile);
            std::ifstream poly("output_polygon.txt");
            std::getline(poly, line); // skip first line
            int i = 0;
            char *temp;
            int pos;
            int x;
            int y;
            std::string sub1;
            std::string sub2;
            while (i < file_name_counter)
            { // extracting the points of the polygon from my output and making the polygon we will use
                std::getline(poly, line);
                pos = line.find(" ");
                sub1 = line.substr(0, pos);
                sub2 = line.substr(pos + 1);
                x = stoi(sub1);
                y = stoi(sub2);
                p.push_back(Point_2(x, y));
                i++;
            }
            get_points(file_name_counter);
            create_chain(file_name_counter);
            Points temp_p = points;
            segments temp_c = chain;
            Polygon_2 temp_po = p;
            for (int k = 0; k < 5; k++)
            {
                L = 1;
                threshold = 0.1;
                std::ofstream outf;
                outf.open("output.txt", std::ofstream::out | std::ofstream::trunc);
                auto start1 = high_resolution_clock::now();
                double ratio = local_search(min_m, outf); // min
                auto stop1 = std::chrono::high_resolution_clock::now();
                auto duration1 = duration_cast<std::chrono::milliseconds>(stop1 - start1); // ton xrono na ton upologizw mesa sthn local kai na ths dinw orisma to cut off an to xepernaei enw trexei na kanei ena flag 1
                if (duration1.count() > cut_off)
                {
                    if (min_m == 1)
                    {
                        ratio == 1;
                    }
                    else
                    {
                        ratio = 0;
                    }
                }
                ratio_array[k] = ratio; // tha exw ola ta ratio
                std::cout << ratio << std::endl;
                points.clear();
                chain.clear();
                p.clear();
                points = temp_p;
                chain = temp_c;
                p = temp_po;
                outf.close();
                if (k == 4 && min_m == 1)
                {
                    double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                    project_out << mo << "\t||";
                    min=std::min_element(ratio_array,ratio_array+5);
                    k = 0;
                    min_m = 2; // max
                }
            }
            double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
            project_out << mo << "\t||";
            project_out <<*min<<"\t||";
            project_out<<*(std::max_element(ratio_array,ratio_array+5))<<"\t||" <<std::endl <<"\t||";
            p.clear();
            points.clear();
            chain.clear();
            poly.close();
            outfile.close();
        }
        for (int j = 1; j <= 3; j++) // incremental random min max 1a 2a 1b 2b local serach min max
        {
            outfile.open("output_polygon.txt", std::ofstream::out | std::ofstream::trunc); // output file
            flagedge = j;                                                                  // random ,min ,max
            int min_m = 1;                                                                 // min
            for (int kl = 1; kl <= 4; kl++)
            {
                flaginit = kl;
                incremental_fun(result2, outfile);
                std::ifstream poly("output_polygon.txt");
                std::getline(poly, line); // skip first line
                int i = 0;
                char *temp;
                int pos;
                int x;
                int y;
                std::string sub1;
                std::string sub2;
                while (i < file_name_counter)
                { // extracting the points of the polygon from my output and making the polygon we will use
                    std::getline(poly, line);
                    pos = line.find(" ");
                    sub1 = line.substr(0, pos);
                    sub2 = line.substr(pos + 1);
                    x = stoi(sub1);
                    y = stoi(sub2);
                    p.push_back(Point_2(x, y));
                    i++;
                }
                get_points(file_name_counter);
                create_chain(file_name_counter);
                Points temp_p = points;
                segments temp_c = chain;
                Polygon_2 temp_po = p;
                for (int k = 0; k < 5; k++)
                {
                    L = 4;
                    threshold = 0.1;
                    std::ofstream outf;
                    outf.open("output.txt", std::ofstream::out | std::ofstream::trunc);
                    auto start1 = high_resolution_clock::now();
                    double ratio = local_search(min_m, outf); // min
                    auto stop1 = std::chrono::high_resolution_clock::now();
                    auto duration1 = duration_cast<std::chrono::milliseconds>(stop1 - start1); // ton xrono na ton upologizw mesa sthn local kai na ths dinw orisma to cut off an to xepernaei enw trexei na kanei ena flag 1
                    if (duration1.count() > cut_off)
                    {
                        if (min_m == 1)
                        {
                            ratio == 1;
                        }
                        else
                        {
                            ratio = 0;
                        }
                    }
                    ratio_array[k] = ratio; // tha exw ola ta ratio
                    std::cout << ratio << std::endl;
                    points.clear();
                    chain.clear();
                    p.clear();
                    points = temp_p;
                    chain = temp_c;
                    p = temp_po;
                    outf.close();
                    if (k == 4 && min_m == 1)
                    {
                        double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                        project_out << mo << "\t||";
                        min=std::min_element(ratio_array,ratio_array+5);
                        k = 0;
                        min_m = 2; // max
                    }
                }
                double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                project_out << mo << "\t||";
                project_out <<*min<<"\t||";
                project_out<<*(std::max_element(ratio_array,ratio_array+5))<<"\t||" <<std::endl <<"\t||";
                p.clear();
                points.clear();
                chain.clear();
                poly.close();
                outfile.close();
            }
        }

        for (int j = 1; j <= 3; j++) // convex hull random min max simulated annealing global
        {
            outfile.open("output_polygon.txt", std::ofstream::out | std::ofstream::trunc); // output file
            flagedge = j;                                                                  // random ,min ,max
            flag_min_max = 1;                                                              // min
            convex_hull_fun(result2, outfile);
            std::ifstream poly("output_polygon.txt");
            std::getline(poly, line); // skip first line
            int i = 0;
            char *temp;
            int pos;
            int x;
            int y;
            std::string sub1;
            std::string sub2;
            while (i < file_name_counter)
            { // extracting the points of the polygon from my output and making the polygon we will use
                std::getline(poly, line);
                pos = line.find(" ");
                sub1 = line.substr(0, pos);
                sub2 = line.substr(pos + 1);
                x = stoi(sub1);
                y = stoi(sub2);
                p.push_back(Point_2(x, y));
                i++;
            }
            get_points(file_name_counter);
            create_chain(file_name_counter);
            Points temp_p = points;
            segments temp_c = chain;
            Polygon_2 temp_po = p;
            for (int k = 0; k < 5; k++)
            {
                L = 500;
                option = 1; // local
                std::ofstream outf;
                outf.open("output.txt", std::ofstream::out | std::ofstream::trunc);
                auto start1 = high_resolution_clock::now();
                double ratio = simulated_annealing(file_name_counter, outf); // min
                auto stop1 = std::chrono::high_resolution_clock::now();
                auto duration1 = duration_cast<std::chrono::milliseconds>(stop1 - start1); // ton xrono na ton upologizw mesa sthn local kai na ths dinw orisma to cut off an to xepernaei enw trexei na kanei ena flag 1
                if (duration1.count() > cut_off)
                {
                    if (flag_min_max == 1)
                    {
                        ratio == 1;
                    }
                    else
                    {
                        ratio = 0;
                    }
                }
                ratio_array[k] = ratio; // tha exw ola ta ratio
                std::cout << ratio << std::endl;
                points.clear();
                chain.clear();
                p.clear();
                points = temp_p;
                chain = temp_c;
                p = temp_po;
                outf.close();
                if (k == 4 && flag_min_max == 1)
                {
                    double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                    project_out << mo << "\t||";
                    min=std::min_element(ratio_array,ratio_array+5);
                    k = 0;
                    flag_min_max = 2; // max
                }
            }
            double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
            project_out << mo << "\t||";
            project_out <<*min<<"\t||";
            project_out<<*(std::max_element(ratio_array,ratio_array+5))<<"\t||" <<std::endl <<"\t||";
            p.clear();
            points.clear();
            chain.clear();
            poly.close();
            outfile.close();
        }

        for (int j = 1; j <= 3; j++) // convex_hull random min max simulated_annealing local min max
        {
            outfile.open("output_polygon.txt", std::ofstream::out | std::ofstream::trunc); // output file
            flagedge = j;                                                                  // random ,min ,max
            flag_min_max = 1;                                                              // min
            convex_hull_fun(result2, outfile);
            std::ifstream poly("output_polygon.txt");
            std::getline(poly, line); // skip first line
            int i = 0;
            char *temp;
            int pos;
            int x;
            int y;
            std::string sub1;
            std::string sub2;
            while (i < file_name_counter)
            { // extracting the points of the polygon from my output and making the polygon we will use
                std::getline(poly, line);
                pos = line.find(" ");
                sub1 = line.substr(0, pos);
                sub2 = line.substr(pos + 1);
                x = stoi(sub1);
                y = stoi(sub2);
                p.push_back(Point_2(x, y));
                i++;
            }
            get_points(file_name_counter);
            create_chain(file_name_counter);
            Points temp_p = points;
            segments temp_c = chain;
            Polygon_2 temp_po = p;
            for (int k = 0; k < 5; k++)
            {
                L = 500;
                option = 2; // local
                std::ofstream outf;
                outf.open("output.txt", std::ofstream::out | std::ofstream::trunc);
                auto start1 = high_resolution_clock::now();
                double ratio = simulated_annealing(file_name_counter, outf); // min
                auto stop1 = std::chrono::high_resolution_clock::now();
                auto duration1 = duration_cast<std::chrono::milliseconds>(stop1 - start1); // ton xrono na ton upologizw mesa sthn local kai na ths dinw orisma to cut off an to xepernaei enw trexei na kanei ena flag 1
                if (duration1.count() > cut_off)
                {
                    if (flag_min_max == 1)
                    {
                        ratio == 1;
                    }
                    else
                    {
                        ratio = 0;
                    }
                }
                ratio_array[k] = ratio; // tha exw ola ta ratio
                std::cout << ratio << std::endl;
                points.clear();
                chain.clear();
                p.clear();
                points = temp_p;
                chain = temp_c;
                p = temp_po;
                outf.close();
                if (k == 4 && flag_min_max == 1)
                {
                    double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                    project_out << mo << "\t||";
                    min=std::min_element(ratio_array,ratio_array+5);
                    k = 0;
                    flag_min_max = 2; // max
                }
            }
            double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
            project_out << mo << "\t||";
            project_out <<*min<<"\t||";
            project_out<<*(std::max_element(ratio_array,ratio_array+5))<<"\t||" <<std::endl <<"\t||";
            p.clear();
            points.clear();
            chain.clear();
            poly.close();
            outfile.close();
        }
        for (int j = 1; j <= 3; j++) // incremental random min max 1a 2a 1b 2b simulated annealing global min max
        {
            outfile.open("output_polygon.txt", std::ofstream::out | std::ofstream::trunc); // output file
            flagedge = j;                                                                  // random ,min ,max
            flag_min_max = 1;                                                              // min
            for (int kl = 1; kl <= 4; kl++)
            {
                flaginit = kl;
                incremental_fun(result2, outfile);
                std::ifstream poly("output_polygon.txt");
                std::getline(poly, line); // skip first line
                int i = 0;
                char *temp;
                int pos;
                int x;
                int y;
                std::string sub1;
                std::string sub2;
                while (i < file_name_counter)
                { // extracting the points of the polygon from my output and making the polygon we will use
                    std::getline(poly, line);
                    pos = line.find(" ");
                    sub1 = line.substr(0, pos);
                    sub2 = line.substr(pos + 1);
                    x = stoi(sub1);
                    y = stoi(sub2);
                    p.push_back(Point_2(x, y));
                    i++;
                }
                get_points(file_name_counter);
                create_chain(file_name_counter);
                Points temp_p = points;
                segments temp_c = chain;
                Polygon_2 temp_po = p;
                for (int k = 0; k < 5; k++)
                {
                    L = 500;
                    option = 1;
                    std::ofstream outf;
                    outf.open("output.txt", std::ofstream::out | std::ofstream::trunc);
                    auto start1 = high_resolution_clock::now();
                    double ratio = simulated_annealing(file_name_counter, outf); // min
                    auto stop1 = std::chrono::high_resolution_clock::now();
                    auto duration1 = duration_cast<std::chrono::milliseconds>(stop1 - start1); // ton xrono na ton upologizw mesa sthn local kai na ths dinw orisma to cut off an to xepernaei enw trexei na kanei ena flag 1
                    if (duration1.count() > cut_off)
                    {
                        if (flag_min_max == 1)
                        {
                            ratio == 1;
                        }
                        else
                        {
                            ratio = 0;
                        }
                    }
                    ratio_array[k] = ratio; // tha exw ola ta ratio
                    std::cout << ratio << std::endl;
                    points.clear();
                    chain.clear();
                    p.clear();
                    points = temp_p;
                    chain = temp_c;
                    p = temp_po;
                    outf.close();
                    if (k == 4 && flag_min_max == 1)
                    {
                        double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                        project_out << mo << "\t||";
                        min=std::min_element(ratio_array,ratio_array+5);
                        k = 0;
                        flag_min_max = 2; // max
                    }
                }
                double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                project_out << mo << "\t||";
                project_out <<*min<<"\t||";
                project_out<<*(std::max_element(ratio_array,ratio_array+5))<<"\t||" <<std::endl <<"\t||";
                p.clear();
                points.clear();
                chain.clear();
                poly.close();
                outfile.close();
            }
        }

        for (int j = 1; j <= 3; j++) // incremental random min max 1a 1b 2a 2b local min max
        {
            outfile.open("output_polygon.txt", std::ofstream::out | std::ofstream::trunc); // output file
            flagedge = j;                                                                  // random ,min ,max
            flag_min_max = 1;                                                              // min
            for (int kl = 1; kl <= 4; kl++)
            {
                flaginit = kl;
                incremental_fun(result2, outfile);
                std::ifstream poly("output_polygon.txt");
                std::getline(poly, line); // skip first line
                int i = 0;
                char *temp;
                int pos;
                int x;
                int y;
                std::string sub1;
                std::string sub2;
                while (i < file_name_counter)
                { // extracting the points of the polygon from my output and making the polygon we will use
                    std::getline(poly, line);
                    pos = line.find(" ");
                    sub1 = line.substr(0, pos);
                    sub2 = line.substr(pos + 1);
                    x = stoi(sub1);
                    y = stoi(sub2);
                    p.push_back(Point_2(x, y));
                    i++;
                }
                get_points(file_name_counter);
                create_chain(file_name_counter);
                Points temp_p = points;
                segments temp_c = chain;
                Polygon_2 temp_po = p;
                for (int k = 0; k < 5; k++)
                {
                    L = 500;
                    option = 2;
                    std::ofstream outf;
                    outf.open("output.txt", std::ofstream::out | std::ofstream::trunc);
                    auto start1 = high_resolution_clock::now();
                    double ratio = simulated_annealing(file_name_counter, outf); // min
                    auto stop1 = std::chrono::high_resolution_clock::now();
                    auto duration1 = duration_cast<std::chrono::milliseconds>(stop1 - start1); // ton xrono na ton upologizw mesa sthn local kai na ths dinw orisma to cut off an to xepernaei enw trexei na kanei ena flag 1
                    if (duration1.count() > cut_off)
                    {
                        if (flag_min_max == 1)
                        {
                            ratio == 1;
                        }
                        else
                        {
                            ratio = 0;
                        }
                    }
                    ratio_array[k] = ratio; // tha exw ola ta ratio
                    std::cout << ratio << std::endl;
                    points.clear();
                    chain.clear();
                    p.clear();
                    points = temp_p;
                    chain = temp_c;
                    p = temp_po;
                    outf.close();
                    if (k == 4 && flag_min_max == 1)
                    {
                        double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                        project_out << mo << "\t||";
                        min=std::min_element(ratio_array,ratio_array+5);
                        k = 0;
                        flag_min_max = 2; // max
                    }
                }
                double mo = (ratio_array[0] + ratio_array[1] + ratio_array[2] + ratio_array[3] + ratio_array[4]) / 5;
                project_out << mo << "\t||";
                project_out <<*min<<"\t||";
                project_out<<*(std::max_element(ratio_array,ratio_array+5))<<"\t||" <<std::endl <<"\t||";
                p.clear();
                points.clear();
                chain.clear();
                poly.close();
                outfile.close();
            }
        }

        if (file_name_counter < 100)
        {
            file_name_counter = file_name_counter + 10;
        }
        else if (file_name_counter < 1000)
        {
            file_name_counter = file_name_counter + 100;
        }
        else if (file_name_counter < 10000)
        {
            file_name_counter = file_name_counter + 1000;
        }
        else if (file_name_counter < 100000)
        {
            file_name_counter = file_name_counter + 10000;
        }
        in_f.close();
    }
}
