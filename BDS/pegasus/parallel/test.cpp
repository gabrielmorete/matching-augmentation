#include <stdio>
#include <omp.h>

using namespace std;

#define NUM_THREADS 6

int cnt = 0;
bool read(string &s, int &my_cnt){
	bool ok = 1;
	#pragma omp critical
	{
		ok = (bool)(cin>>s);
		if (ok)
			cnt++;
		my_cnt = cnt;
	}
	return ok;
}

signed main(){
	string s;

    #pragma omp parallel num_threads(NUM_THREADS)\
    	private(s)\
    	shared(cnt)
    {
    	int my_cnt;	
		while (read(s, my_cnt)){
			int n = omp_get_thread_num(); 
			#pragma omp critical
			{
				cout<<n<<' '<<s<<' '<<my_cnt<<endl;
			}
		}
	}	

}