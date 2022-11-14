#include <iterator>
#include <map>
#include <cstdio>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
//#include <omp.h>
//#include <string>
#include <unordered_map>


using namespace std;




struct event{
	
	double time;
	long id1;
	long id2;
	double m1;
	double m2;
	double ecc;
	double a;


};
	
typedef unordered_multimap<long long, event> hmap;
typedef long long key_type;
const auto relative_difference_factor = 0.0005;  // 0.05%
hmap pos_map;
hmap neg_map;

void print_event(event m)
{
	printf("Time is %8.8f, id1= %d, id2= %d, m1= %4.8f, m2= %4.8f, ecc= %4.4f, a= %4.4f\n", m.time,m.id1,m.id2,m.m1,m.m2,m.ecc,m.a);
}


key_type hash_func(long id1,long id2)
{	key_type hs;
//string hs;
//hs=to_string(id1) + "," + to_string(id2);
	hs=id2+((id1 + id2 +1)*(id1+id2)/2);
	return hs;

};

bool check_key(hmap m, key_type key)
{	//printf("Number of keys: %d\n",m.count(key));
	if (m.count(key)==0)
		return false;
	return true;

}

void save_csv(hmap map,const char* storage_path)
{	FILE * outfile;
	outfile=fopen(storage_path,"w");
	for (auto it = map.begin(); it != map.end(); it++) 
		{
         event f=it->second;
		 fprintf(outfile,"%8.3f,%7d,%7d,%7.12f,%7.12f,%8.10f,%8.10f\n",f.time, f.id1, f.id2, f.m1, f.m2, f.a,f.ecc);
		}
	fclose(outfile);
}



hmap data_reduction()
{
	event min_pos_ev;
	bool min_found;
	float min_time_diff;

	hmap matching_map;
//	#pragma omp for
	for (hmap::iterator it=neg_map.begin(); it!=neg_map.end(); it++) 
//	for(int i =0; i<neg_map.size(); i++)
	{	//printf("Negative event\n");
//		auto it=neg_map.begin();		
//		advance(it,i);
		const key_type key=it->first;
		
		event neg_ev=it->second;
		//print_event(neg_ev);
		if (check_key(pos_map, key)) 
		{	//printf("Positive events\n");
			auto pos_itr=pos_map.equal_range(key);	
			min_found=false;
			min_time_diff=neg_ev.time;
			for (auto pos_it=pos_itr.first;pos_it!=pos_itr.second;pos_it++)
			{	
				event pos_ev=pos_it->second;
				if (pos_ev.time<=neg_ev.time && min_time_diff>neg_ev.time-pos_ev.time)
				{
					min_time_diff=neg_ev.time-pos_ev.time;
					min_pos_ev=pos_ev;
					min_found=true;
			//		print_event(pos_ev);
				}
			}
			const auto greater_magnitude = max(abs(min_pos_ev.a),abs(neg_ev.a));
			if (min_found && abs(min_pos_ev.a-neg_ev.a)<relative_difference_factor*greater_magnitude) 
			{
				matching_map.insert({key,neg_ev});
			}
		}
	//	printf("End of neg event\n");
	}
	return matching_map;

}


void read_csv()
{

	FILE * orb_file;
	orb_file=fopen("orb_elements.csv","r");
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
		key_type key;
		key=hash_func(temp_event.id1,temp_event.id2);
		if (*split[0]=='+') 
		{			
			pos_map.insert({key,temp_event});
		}
		else 
		{
			neg_map.insert({key,temp_event});
		}

	};
	fclose(orb_file);
	
}




int main(){

	read_csv();
	save_csv(pos_map, "test_pos.csv");
	save_csv(neg_map, "test_neg.csv");
	hmap matching_map=data_reduction();
	printf("The load factor of the new dataset is: %2.2f",matching_map.load_factor());
	save_csv(matching_map, "test_reduced.csv");







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
