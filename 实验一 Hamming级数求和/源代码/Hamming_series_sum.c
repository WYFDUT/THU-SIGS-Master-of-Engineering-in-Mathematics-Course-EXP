#include<stdio.h>
#include<math.h>
#include<time.h>

void HammingSeriesSumExpansionSolver(double res[]){
    double x = 0.0, sum = 0.0;
    int k = 1, index = 0;
    while(x <= 300.0){
        sum = 0.0;
        k = 1;
        while(k <= 60000){
            sum += 1.0/(((double)k)*((double)k+1.0)*((double)k+2.0)*((double)k+x));
            k += 1;
        }
        sum = 1.0 + (1.0-x) * (2.0-x) * sum + (1.0-x)/4.0;
        //printf("%3.1f  %2.11f\n", x, sum);
        res[index] = sum;
        x += 0.1;
        index += 1;
    }
}

void HammingSeriesSumRecursiveSolver(double res[]){
    int i, k;
	double x[3001];
	for (i = 0; i < 3001; i++)			// set 0.0 ~ 300.0 in x[ ]
		x[i] = (double)i * 0.10;

	for (i = 0; i < 10; i++) {			// calculate 0.0 ~ 0.9
		res[i] = 0;
		for (k = 1; k < 2000; k++) {	// up bound is got from testing
			res[i] += 1.0 / ((double)k * (double)(k+1) * (double)(k+2) * ((double)k + x[i]));
		}
		res[i] = (res[i] * (2.0-x[i]) + 0.25) * (1.0 - x[i]) + 1.0;
	}

	for (i = 10; i <= 3000; i++) {		// calculate 1.0 ~ 299.0
		res[i] = ( 1.0/x[i] + (x[i]-1.0) * res[i-10] ) / x[i];
	}
}


int main(){
    double sum[3001];
    int begintime, endtime;
    FILE *fp;
    begintime=clock();	//计时开始
    HammingSeriesSumRecursiveSolver(sum);
    endtime = clock();	//计时结束
    fp = fopen("result.txt", "w");
    for(int i = 0; i<3001; i++){
        printf("%3.1f  %2.10f\n", (double)i/10, sum[i]);
        fprintf(fp, "%3.1f  %2.10f\n", (double)i/10, sum[i]);
    }
	printf("\n\nRunning Time: %dms\n", endtime-begintime);
    fclose(fp);
    system("pause");
    return 0;
}