// Code to analyse the orb_elements file to find interactions of binaries in the system, mainly catches hyperbolic interaction.


#include <iterator>
#include <map>
#include <cstdio>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <omp.h>
#include <string>
//#include <mpi.h>
#include <algorithm>
#include <string>
#include <vector>
#include <filesystem>

using namespace std;


const int NUM_THREADS=10;

struct event{
	
	double time;
	long id1;
	long id2;
	double m1;
	double m2;
	double ecc;
	double a;


};
	
const auto relative_difference_factor = 0.0005;  // 0.05%
void print_event(event m)
{
	printf("Time is %8.8f, id1= %d, id2= %d, m1= %4.8f, m2= %4.8f, ecc= %4.4f, a= %4.4f\n", m.time,m.id1,m.id2,m.m1,m.m2,m.ecc,m.a);
}

void save_csv(vector<event> map,const char* storage_path, int file_num)
{	FILE * outfile;
	if (file_num!=0){
	outfile=fopen(storage_path,"a");
	}
	else {
	outfile=fopen(storage_path,"w");
	}
	for (auto it = map.begin(); it != map.end(); it++) 
		{
         event f=*it;
		 fprintf(outfile,"%8.3f,%7d,%7d,%7.12f,%7.12f,%8.10f,%8.10f\n",f.time, f.id1, f.id2, f.m1, f.m2, f.a,f.ecc);
		}
	fclose(outfile);
}
vector <vector <event>::iterator> find_index( vector <event> map,long fid1,long fid2)
{	vector <event>::iterator iter = begin(map);
	vector <vector<event>::iterator> iterator_arr;
	while((iter = find_if(iter, map.end(), [=](event found_event){
if(found_event.id1 == fid1 &&found_event.id2==fid2)
	{
		return true;
	}
else{
	return false;
	}
}))!=end(map))
			{
				if(iter!=map.end()){
				iterator_arr.push_back(iter);
				}
				iter++;
			}


	return iterator_arr;
}


vector<event> data_reduction(vector <event> &pos_map, vector <event> &neg_map)
{
	vector<event> matching_map[NUM_THREADS];

	int len_map=neg_map.size();
	int n_per_thread = len_map/NUM_THREADS;
	int i;

#pragma omp parallel num_threads(NUM_THREADS) shared(matching_map, pos_map,neg_map) private(i) 
	{	
		int count=0;

	//#pragma omp for schedule(static,1)
	for(i=omp_get_thread_num()*n_per_thread; i<(omp_get_thread_num()+1)*n_per_thread; i++)
	{	
		//printf("Negative event\n");
	//	if (omp_get_thread_num()==0)
	//	printf("%d\n",i);
		
		auto it=neg_map.begin();		
		advance(it,i);
		const event neg_ev=*it;
		const long fid1=neg_ev.id1;
		const long fid2=neg_ev.id2;
		event min_pos_ev;
		bool min_found;
		float min_time_diff;
		auto pos_itr_arr=find_index(pos_map,fid1,fid2);
		//printf("size of pos_itr_arr is %d\n",pos_itr_arr.size());
		//print_event(neg_ev);

		//printf("index %d\n",i);
		if (pos_itr_arr.size()!=0) 
		{	//printf("Positive events\n");
			//auto pos_itr=pos_map.equal_range(key);	
			count++;
			min_found=false;
			min_time_diff=neg_ev.time;
			for (auto pos_it_point=pos_itr_arr.begin(); pos_it_point<pos_itr_arr.end();pos_it_point++)
			{	auto pos_it=*pos_it_point;
				if(pos_it<pos_map.end())
				{	
				event pos_ev=*pos_it;

					if (pos_ev.time<=neg_ev.time && min_time_diff>neg_ev.time-pos_ev.time)
					{
					min_time_diff=neg_ev.time-pos_ev.time;
					min_pos_ev=pos_ev;
					min_found=true;

					}
	
				}


			}
			if (min_found && abs(min_pos_ev.a-neg_ev.a)<relative_difference_factor* max(abs(min_pos_ev.a),abs(neg_ev.a))) 
			{
				{
				matching_map[omp_get_thread_num()].push_back(neg_ev);
				}
			}
		}

	}
	printf("index for thread %d is %d-%d\n",omp_get_thread_num(),omp_get_thread_num()*n_per_thread,(omp_get_thread_num()+1)*n_per_thread);
		printf("count of thread %d is %d\n",omp_get_thread_num(),count);
	printf("size of matching map for thread %d is %d\n",omp_get_thread_num(),matching_map[omp_get_thread_num()].size());
	
}

//printf("size of matching map total is %d\n",int(matching_map[0].size())+int(matching_map[1].size()));

vector <event> flattened_matching_map;

for(int th=0; th<NUM_THREADS;th++)
{
	for(auto len=matching_map[th].begin();len<matching_map[th].end();len++)
	{	
		flattened_matching_map.push_back(*len);
	}

}
return flattened_matching_map;

}


void read_csv(vector <event> & pos_map, vector <event> &neg_map,int file_num)
{

	FILE * orb_file;
	string file_num_append;
	file_num_append=to_string(file_num);
	string file_name="split/orb_elements_"+file_num_append+".csv";
	const char* file_name_char = file_name.c_str();
	orb_file=fopen(file_name_char,"r");
	char content[1024];
	while(fgets(content,1024,orb_file))
	{
		char *split[8];
		int num_token=0;
		//printf("%s\n",content);
		char *token = strtok(content, ",");
		while(token!=NULL)
		{
			split[num_token++]= token;
			token=strtok(NULL,",");
		};
	
		struct event temp_event;
		temp_event.time=stod(split[1]);
		temp_event.id1=stol(split[2]);
		temp_event.id2=stol(split[3]);
		temp_event.m1=stod(split[4]);
		temp_event.m2=stod(split[5]);
		temp_event.ecc=stod(split[6]);
		temp_event.a=stod(split[7]);
		if (*split[0]=='+') 
		{			
			pos_map.push_back(temp_event);
		}
		else 
		{
			neg_map.push_back(temp_event);
		}

	};
	fclose(orb_file);
	
}




int main(){
int file_count=0;

for (auto & entry : std::filesystem::directory_iterator("split/"))
        file_count++;
printf("number of files is %d\n",file_count);

for(int file_num=0;file_num<file_count;file_num++){
	vector<event> pos_map;
	vector<event> neg_map;
	read_csv(pos_map,neg_map,file_num);
	save_csv(pos_map, "pos.csv",file_num);
	save_csv(neg_map, "neg.csv",file_num);

	auto matching_map=data_reduction(pos_map,neg_map);
	save_csv(matching_map, "reduced.csv",file_num);
	printf("file num %d completed\n",file_num);
}

/*	key_type testkey= hash_func(41325,  41326);
	event f;
	if (check_key(pos_map, testkey))
	{   auto itr = pos_map.equal_range(testkey);
		for (auto it = itr.first; it != itr.second; it++) 
		{
         f=it->second;
		 print_event(f);
		}	

	}
	else {
	printf("key not present");
	}

	//f=neg_map.at(testkey);
	//printf("%2.10f",f.m1);
*/


	return 0;
}
