#define _CRT_SECURE_NO_WARNINGS     // sprintf is Err in VC++

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>
#include <ctype.h>



#include <cuda_runtime.h>

#include <cstring>
#include <cstdlib>
#include <vector>

#include <string>
#include <iostream>
#include <stdio.h>

#include "caffe/caffe.hpp"
#include "caffe/util/io.hpp"
#include "caffe/blob.hpp"

using namespace caffe;
using namespace std;



const int BLACK = 1;
const int WHITE = 2;


int BubSort(int x[], int n);

double komi = 6.5;

int use_time_control = 1;    // 0 ... no time limit

#define B_SIZE     19
#define WIDTH      (B_SIZE + 2)
#define BOARD_MAX  (WIDTH * WIDTH)



#if defined(_MSC_VER)
typedef unsigned __int64 uint64;
#define PRIx64  "I64x"
#else
#include <stdint.h>
#include <sys/time.h>
#include <unistd.h>
typedef uint64_t uint64;  // Linux
#define PRIx64  "llx"
#endif

#define HASH_KINDS 4      // 1...black, 2...white, 3...ko
#define HASH_KO 3
uint64 hashboard[BOARD_MAX][HASH_KINDS];
uint64 hashcode = 0;




int dir4[4] = { +1,+WIDTH,-1,-WIDTH };
int dir8[8] = { +1,+WIDTH,-1,-WIDTH,  +1 + WIDTH, +WIDTH - 1, -1 - WIDTH, -WIDTH + 1 };

int ko_z;

#define MAX_MOVES 9000
int record[MAX_MOVES];
double record_time[MAX_MOVES];
int moves = 0;

int all_playouts = 0;
int flag_test_playout = 0;

#define D_MAX 9000
int path2[D_MAX];
int depth;

int board_area_sum[BOARD_MAX];
int board_winner[2][BOARD_MAX];
int winner_count[2];


int criticality[BOARD_MAX];


void prt(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	//{ FILE *fp = fopen("out.txt","a"); if ( fp ) { vfprt( fp, fmt, ap ); fclose(fp); } }
	vfprintf(stderr, fmt, ap);
	va_end(ap);
}
void send_gtp(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfprintf(stdout, fmt, ap);
	va_end(ap);
}



unsigned long rand_xorshift128() {  // 2^128-1 
	static unsigned long x = 123456789, y = 362436069, z = 521288629, w = 88675123;
	unsigned long t;
	t = (x ^ (x << 11)) & 0xffffffff;
	x = y; y = z; z = w; return(w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

uint64 rand64()
{
	unsigned long r1 = rand_xorshift128();
	unsigned long r2 = rand_xorshift128();
	uint64 r = ((uint64)r1 << 32) | r2;
	return r;
}

void prt_code64(uint64 r)
{
	//prt("%016" PRIx64,r);
	prt("%08x%08x", (int)(r >> 32), (int)r);
};

void make_hashboard()
{
	int z, i;
	for (z = 0; z<BOARD_MAX; z++) {
		//  prt("[%3d]=",z); 
		for (i = 0; i<HASH_KINDS; i++) {
			hashboard[z][i] = rand64();
			//    prt_code64(hashboard[z][i]); prt(",");
		}
		//  prt("\n");
	}
}

void hash_pass()
{
	hashcode = ~hashcode;
}
void hash_xor(int z, int color)
{
	hashcode ^= hashboard[z][color];
}



int get_z(int x, int y)
{
	return y*WIDTH + x;  // 1<= x <=9, 1<= y <=9
}

int get81(int z)            // for display only
{
	int y = z / WIDTH;
	int x = z - y*WIDTH;    // 106 = 9*11 + 7 = (x,y)=(7,9) -> 79
	if (z == 0) return 0;
	return x * 100 + y;        // x*100+y for 19x19
}

// don't call twice in same sentence. like prt("z0=%s,z1=%s\n",get_char_z(z0),get_char_z(z1));
char *get_char_z(int z)
{
	int x, y, ax;
	static char buf[16];
	sprintf(buf, "pass");
	if (z == 0) return buf;
	y = z / WIDTH;
	x = z - y*WIDTH;
	ax = x - 1 + 'A';
	if (ax >= 'I') ax++;  // from 'A' to 'T', excluding 'I'
	sprintf(buf, "%c%d", ax, B_SIZE + 1 - y);
	return buf;
}

int flip_color(int col)
{
	return 3 - col;
}

/*
int board[BOARD_MAX] = {
  3,3,3,3,3,3,3,3,3,3,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,0,0,0,0,0,0,0,0,0,3,
  3,3,3,3,3,3,3,3,3,3,3
};


int board[BOARD_MAX] = {
	3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,

};
*/

int board[BOARD_MAX] = {
	3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
	3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,

};



int board2[19*19] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,2,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,1,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,1,1,2,1,0,0,0,0,0,0,0,0,0,1,0,2,0,0,
	0,1,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

float board3[19*19] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};


//const int B_SIZE = 19;

int xy2pos(int x,int y)
{
	return y*B_SIZE+x;
}


Net<float> *p_caffe_net;

void initCNN()
{
	//set_phase(TEST);
	Caffe::set_mode(Caffe::CPU);
	//Caffe::set_mode(Caffe::GPU);

	string net_src = "movepredict.prototxt";
//	Net<float> caffe_net(net_src, caffe::TEST);  //get the net

//	p_caffe_net = &caffe_net;
	p_caffe_net = new Net<float>(net_src, caffe::TEST);
   
	string traied_net = "movepredict.caffemodel";
	p_caffe_net->CopyTrainedLayersFrom(traied_net);
}
void testCNN(int col2,float *q)
{
	int col = col2;
	const int size = B_SIZE;
	float result[size*size];
	float *data = new float[2*size*size];
	//fprintf(stderr,"2\n");
	if (col==BLACK) {
		for (int j=0;j<size;j++) for (int k=0;k<size;k++) {
//			board2[xy2pos(j,k)]=board[get_z(j + 1, k + 1)];	
//fprintf(stderr,"%d %d %d\n",i,j,k);
			if (board2[xy2pos(j,k)]==BLACK) {
//			if (board[get_z(j + 1, k + 1)]==BLACK) {
				data[j*size+k]=1.0;
				data[size*size+size*j+k]=0.0;
			} else if (board[get_z(j + 1, k + 1)]==WHITE) {
				data[j*size+k]=0.0;
				data[size*size+size*j+k]=1.0;
			} else {
				data[j*size+k]=0.0;
				data[size*size+size*j+k]=0.0;
			}
		}
//		++p;
	}
	if (col==WHITE) {
		for (int j=0;j<size;j++) for (int k=0;k<size;k++) {
//			board2[xy2pos(j,k)]=*p;	

			//fprintf(stderr,"%d %d %d\n",i,j,k);
//			board2[xy2pos(j,k)]=board[get_z(j + 1, k + 1)];	
			if (board2[xy2pos(j,k)]==BLACK) {
//			if (board[get_z(j + 1, k + 1)]==BLACK) {
				data[j*size+k]=0.0;
				data[size*size+size*j+k]=1.0;
			} else if (board[get_z(j + 1, k + 1)]==WHITE) {
				data[j*size+k]=1.0;
				data[size*size+size*j+k]=0.0;
			} else {
				data[j*size+k]=0.0;
				data[size*size+size*j+k]=0.0;
			}
		}
//		++p;
	}

	Blob<float> *b=new Blob<float>(1,2,size,size);
	b->set_cpu_data(data);
	vector<Blob<float>*> bottom;
	bottom.push_back(b);
	const vector<Blob<float>*>& rr = p_caffe_net->Forward(bottom);

	float sum = 0;
	for (int j=0;j<B_SIZE;j++) {
		for (int k=0;k<B_SIZE;k++) {
			float a = rr[0]->cpu_data()[j*B_SIZE+k];
			sum += a;
//			fprintf(stderr,"%5.3f ",a);
//		prt("%5.3f ", a);
			*q=(a*1000);//\94z\97\F1\83{\81[\83h3\82Ɋm\97\A6\82\F0\82\A2\82\EA\82\E9
//			board3[xy2pos(j,k)]=(a*1000);//\94z\97\F1\83{\81[\83h3\82Ɋm\97\A6\82\F0\82\A2\82\EA\82\E9
			++q;
		}
//		fprintf(stderr,"\n");
//		prt("\n");
	}
//	fprintf(stderr,"sum=%.3f\n",sum);
//		prt("%.3f\n",sum);

	for (int i=0;i<size*size;i++) {
		result[i]=rr[0]->cpu_data()[i];
		if (result[i]<0.00001) result[i]=0.00001;
	}
	delete[] data;
	delete b;
}


//int main3(int argc, char** argv) {
int polnet(int col2){
//			prt("polnet1");
	initCNN();
//	testCNN();
//			prt("polnet2");

//	double ct1 = clock();
//	for (i=0;i<10;i++) {
		int x,y,tz;

		for (y=0;y<B_SIZE;y++) {
			for (x=0;x<B_SIZE;x++) {
//			int c = rand() % 3;
//			int tz = get_z(x,y);
//			int c = *p//board[tz];
//			board[tz]=*p;
//			board2[xy2pos(x,y)] = board[tz];	
			board3[xy2pos(x,y)] = 0.0;
//		++p;
	}
}

		testCNN(col2,board3);
//	}

//	double t = (clock() - ct1)/CLOCKS_PER_SEC;
//	prt("i=%d,t=%.1f, %f\n",i,t, t/i);

			prt("probability=");
			prt("\n");

		for (y = 0; y<B_SIZE; y++) {
			for (x = 0; x<B_SIZE; x++) {
//			int z = get_z(x + 1, y + 1);
			prt("%5.3f ", board3[xy2pos(x,y)]);
		}
			prt("\n");

	}




	return 0;
}



void print_board()
{
	int x, y;
	const char *str[4] = { ".","X","O","#" };
	int played_z = 0;
	int color = 0;
	if (moves > 0) {
		played_z = record[moves - 1];
		color = board[played_z];
	}

	prt("   ");
	//for (x=0;x<B_SIZE;x++) prt("%d",x+1);
	for (x = 0; x<B_SIZE; x++) prt("%c", 'A' + x + (x>7));
	prt("\n");
	for (y = 0; y<B_SIZE; y++) {
		//  prt("%2d ",y+1);
		prt("%2d ", B_SIZE - y);
		for (x = 0; x<B_SIZE; x++) {
			prt("%s", str[board[get_z(x + 1, y + 1)]]);
		}
		if (y == 4) prt("  ko_z=%s,moves=%d", get_char_z(ko_z), moves);
		if (y == 7) prt("  play_z=%s, color=%d", get_char_z(played_z), color);

		prt("\n");
	}
}

void print_board_area()
{
	int x, y;
	int all = all_playouts;
	if (all == 0) all = 1;

	prt("board_area_sum\n   ");
	for (x = 0; x<B_SIZE; x++) prt("   %c", 'A' + x + (x>7));
	prt("\n");
	for (y = 0; y<B_SIZE; y++) {
		prt("%2d ", B_SIZE - y);
		for (x = 0; x<B_SIZE; x++) {
			int sum = board_area_sum[get_z(x + 1, y + 1)];
			double per = 100.0 * sum / all;
			prt("%4.0f", per);
		}

		prt("\n");
	}
}


double get_criticality(int z)
{
	double all = all_playouts + 1;
	double v = board_winner[0][z] + board_winner[1][z];
	double per = v / all;
	double bp = (double)board_winner[0][z] * winner_count[0] / (all * all);
	double wp = (double)board_winner[1][z] * winner_count[1] / (all * all);
	double criticality = (per - (bp + wp));
	return criticality;
}

void print_criticality()
{
	int x, y;

	prt("criticality\n  ");
	for (x = 0; x<B_SIZE; x++) prt("    %c", 'A' + x + (x>7));
	prt("\n");
	for (y = 0; y<B_SIZE; y++) {
		prt("%2d ", B_SIZE - y);
		for (x = 0; x<B_SIZE; x++) {
			double crt = get_criticality(get_z(x + 1, y + 1));
			prt("%5.2f", crt);
		}
		prt("\n");
	}
}

int check_board[BOARD_MAX];

void count_liberty_sub(int tz, int color, int *p_liberty, int *p_stone)
{
	int z, i;

	check_board[tz] = 1;     // search flag
	(*p_stone)++;            // number of stone
	for (i = 0; i<4; i++) {
		z = tz + dir4[i];
		if (check_board[z]) continue;
		if (board[z] == 0) {
			check_board[z] = 1;
			(*p_liberty)++;      // number of liberty
		}
		if (board[z] == color) count_liberty_sub(z, color, p_liberty, p_stone);
	}
}

void count_liberty(int tz, int *p_liberty, int *p_stone)
{
	int i;
	*p_liberty = *p_stone = 0;
	for (i = 0; i<BOARD_MAX; i++) check_board[i] = 0;
	count_liberty_sub(tz, board[tz], p_liberty, p_stone);
}

void take_stone(int tz, int color)
{
	int z, i;

	hash_xor(tz, color);
	board[tz] = 0;
	for (i = 0; i<4; i++) {
		z = tz + dir4[i];
		if (board[z] == color) take_stone(z, color);
	}
}

const int FILL_EYE_ERR = 1;
const int FILL_EYE_OK = 0;

// put stone. success returns 0. in playout, fill_eye_err = 1.
int put_stone(int tz, int color, int fill_eye_err)
{
	int around[4][3];
	int un_col = flip_color(color);
	int space = 0;
	int wall = 0;
	int mycol_safe = 0;
	int capture_sum = 0;
	int ko_maybe = 0;
	int liberty, stone;
	int i;

	if (tz == 0) {  // pass
		if (ko_z != 0) hash_xor(ko_z, HASH_KO);
		ko_z = 0;
		hash_pass();
		return 0;
	}

	// count 4 neighbor's liberty and stones.
	for (i = 0; i<4; i++) {
		int z, c, liberty, stone;
		around[i][0] = around[i][1] = around[i][2] = 0;
		z = tz + dir4[i];
		c = board[z];  // color
		if (c == 0) space++;
		if (c == 3) wall++;
		if (c == 0 || c == 3) continue;
		count_liberty(z, &liberty, &stone);
		around[i][0] = liberty;
		around[i][1] = stone;
		around[i][2] = c;
		if (c == un_col && liberty == 1) { capture_sum += stone; ko_maybe = z; }
		if (c == color  && liberty >= 2) mycol_safe++;
	}

	if (capture_sum == 0 && space == 0 && mycol_safe == 0) return 1; // suicide
	if (tz == ko_z) return 2; // ko
	if (wall + mycol_safe == 4 && fill_eye_err) return 3; // eye
	if (board[tz] != 0) return 4;

	for (i = 0; i<4; i++) {
		int lib = around[i][0];
		int c = around[i][2];
		if (c == un_col && lib == 1 && board[tz + dir4[i]] != 0) {
			take_stone(tz + dir4[i], un_col);
		}
	}

	board[tz] = color;
	hash_xor(tz, color);
	hash_pass();
	if (ko_z != 0) hash_xor(ko_z, HASH_KO);

	count_liberty(tz, &liberty, &stone);
	if (capture_sum == 1 && stone == 1 && liberty == 1) {
		ko_z = ko_maybe;
		hash_xor(ko_z, HASH_KO);
	}
	else {
		ko_z = 0;
	}
	return 0;
}




int label_min(int a, int b, int c, int d)
{
	int temp_min = 500;

	if (a<temp_min && a>0)temp_min = a;
	if (b<temp_min && b>0)temp_min = b;
	if (c<temp_min && c>0)temp_min = c;
	if (d<temp_min && d>0)temp_min = d;

	return temp_min;
}


// put stone. success returns 0. in playout, fill_eye_err = 1.
int argon_check1(int tz)
{
	int around[4][3];
	//	int un_col = flip_color(color);
	int space = 0;
	int wall = 0;
	int mycol_safe = 0;
	int capture_sum = 0;
	int ko_maybe = 0;
	int liberty, stone;
	int i;

	int label1[WIDTH][WIDTH];
	int label2[WIDTH][WIDTH];
	int number = 1;
	int jinti1[BOARD_MAX];
	int label_num[BOARD_MAX];

	if (board[tz] != 0) return 4;

	if (tz < 1) {  // pass

//		if (ko_z != 0) hash_xor(ko_z, HASH_KO);
//		ko_z = 0;
//		hash_pass();
		return 0;
	}

		if (tz >BOARD_MAX) {  // pass
		return 0;
	}

	// count 4 neighbor's liberty and stones.
	for (i = 0; i<4; i++) {
		int z, c, liberty, stone;
		around[i][0] = around[i][1] = around[i][2] = 0;
		z = tz + dir4[i];
		c = board[z];  // color
		if (c == 0) space++;
		if (c == 3) wall++;
		if (c == 0 || c == 3) continue;
		count_liberty(z, &liberty, &stone);
		around[i][0] = liberty;
		around[i][1] = stone;
		around[i][2] = c;
		//	if (c == un_col && liberty == 1) { capture_sum += stone; ko_maybe = z; }
		//	if (c == color  && liberty >= 2) mycol_safe++;
	}

//	if (capture_sum == 0 && space == 0 && mycol_safe == 0) return 1; // suicide
//	if (tz == ko_z) return 2; // ko
							  //	if (wall + mycol_safe == 4 && fill_eye_err) return 3; // eye
//	if (board[tz] != 0) return 4;



	//argon code 1 \83A\83\8B\83S\83\93\81@\83G\83N\83X\83g\83\89\81@\83R\81[\83h
	for (i = 0; i<BOARD_MAX; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // 106 = 9*11 + 7 = (x,y)=(7,9) -> 79

		label2[x][y] = 0;
	}



	for (i = 0; i<BOARD_MAX; i++) {
		jinti1[i] = 8;

	}


	for (i = 0; i<BOARD_MAX; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // 106 = 9*11 + 7 = (x,y)=(7,9) -> 79
		label1[x][y] = board[i];
		//			        prt("%d ",label1[x][y]);

	}




	for (i = 0; i<BOARD_MAX; i++) {

		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // labeling


		if (label1[x][y] == 0 && label2[x][y] == 0) {
			if (label2[x - 1][y]>0 || label2[x][y - 1]>0 || label2[x + 1][y]>0 || label2[x][y + 1]>0)label2[x][y] = label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]);
			else if (label2[x + 1][y] == 0 && label2[x - 1][y] == 0 && label2[x][y + 1] == 0 && label2[x][y - 1] == 0) {
				label2[x][y] = number;
				number = number + 1;
			}

		}
	}


	for (i = 0; i<BOARD_MAX; i++) {
		int j = 0;
		j = BOARD_MAX -1 - i;
		int y = j / WIDTH;
		int x = j - (j / WIDTH)*WIDTH;    // labeling

		if (label1[x][y] == 0) {
			if (label2[x][y]>label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]) && (label2[x - 1][y]>0 || label2[x][y - 1]>0 || label2[x + 1][y]>0 || label2[x][y + 1]>0))label2[x][y] = label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]);
		}
	}



	for (i = 0; i<BOARD_MAX; i++) {
		int x = i / WIDTH;
		int y = i - (i / WIDTH)*WIDTH;    // labeling

		if (label1[x][y] == 0) {
			if (label2[x][y]>label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]) && (label2[x - 1][y]>0 || label2[x][y - 1]>0 || label2[x + 1][y]>0 || label2[x][y + 1]>0))label2[x][y] = label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]);
		}
	}


	for (i = 0; i<BOARD_MAX; i++) {
		int j = 0;
		j = BOARD_MAX -1 - i;
		int x = j / WIDTH;
		int y = j - (j / WIDTH)*WIDTH;    // labeling

		if (label1[x][y] == 0) {
			if (label2[x][y]>label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]) && (label2[x - 1][y]>0 || label2[x][y - 1]>0 || label2[x + 1][y]>0 || label2[x][y + 1]>0))label2[x][y] = label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]);
		}
	}


	for (i = 0; i<BOARD_MAX; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // labeling

		if (label1[x][y] == 0) {
			if (label2[x][y]>label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]) && (label2[x - 1][y]>0 || label2[x][y - 1]>0 || label2[x + 1][y]>0 || label2[x][y + 1]>0))label2[x][y] = label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]);
		}
	}


/*

	for (i = 0; i<BOARD_MAX; i++) {

		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // labeling


		if (label1[x][y] == 0 && label2[x][y] == 0) {
			if (label2[x - 1][y]>0 || label2[x][y - 1]>0 || label2[x + 1][y]>0 || label2[x][y + 1]>0)label2[x][y] = label_min(label2[x - 1][y], label2[x + 1][y], label2[x][y - 1], label2[x][y + 1]);
			else if (label2[x + 1][y] == 0 && label2[x - 1][y] == 0 && label2[x][y + 1] == 0 && label2[x][y - 1] == 0) {
				label2[x][y] = number;
				number = number + 1;
			}

		}
	}


	for (i = 0; i<BOARD_MAX; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // labeling

		if ( label2[x][y] == 0 && label2[x + 1][y] == 0 && label2[x - 1][y] == 0 && label2[x][y + 1] == 0 && label2[x][y - 1] == 0
			) {
			label2[x][y] = number;
			number = number + 1;
		}
		if ( label2[x][y] == 0 && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] == label2[x - 1][y];
		if ( label2[x][y] == 0 && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] == label2[x][y - 1];
		if ( label2[x][y] == 0 && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] == label2[x + 1][y];
		if ( label2[x][y] == 0 && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] == label2[x][y + 1];

	}


	for (i = BOARD_MAX-1; i <1; i--) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // 

		if (label1[x][y] == 0 && label2[x][y]>label2[x - 1][y] && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] = label2[x - 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y - 1] && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] = label2[x][y - 1];
		if (label1[x][y] == 0 && label2[x][y]>label2[x + 1][y] && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] = label2[x + 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y + 1] && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] = label2[x][y + 1];
	}


	for (i = 0; i<BOARD_MAX; i++) {
		int x = i / WIDTH;
		int y = i - (i / WIDTH)*WIDTH;    // 

		if (label1[x][y] == 0 && label2[x][y]>label2[x - 1][y] && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] = label2[x - 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y - 1] && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] = label2[x][y - 1];
		if (label1[x][y] == 0 && label2[x][y]>label2[x + 1][y] && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] = label2[x + 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y + 1] && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] = label2[x][y + 1];
	}

	for (i = BOARD_MAX-1; i <1; i--) {
	int x = i / WIDTH;
	int y = i - (i / WIDTH)*WIDTH;    //

	if (label1[x][y] == 0 && label2[x][y]>label2[x - 1][y] && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] = label2[x - 1][y];
	if (label1[x][y] == 0 && label2[x][y]>label2[x][y - 1] && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] = label2[x][y - 1];
	if (label1[x][y] == 0 && label2[x][y]>label2[x + 1][y] && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] = label2[x + 1][y];
	if (label1[x][y] == 0 && label2[x][y]>label2[x][y + 1] && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] = label2[x][y + 1];
	}

	for (i = 0; i<BOARD_MAX; i++) {
	int y = i / WIDTH;
	int x = i - (i / WIDTH)*WIDTH;    // labeling

	if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] == label2[x - 1][y];
	if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] == label2[x][y - 1];
	if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] == label2[x + 1][y];
	if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] == label2[x][y + 1];

	}
		*/



	/*
	for (i=0;i<121;i++) {
	int x = i / WIDTH;
	int y = i - (i / WIDTH)*WIDTH;    //
	if (label2[x][y]!=0 && label2[x][y]>label2[x+1][y] && label2[x+1][y]>0)label2[x][y]=label2[x+1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x-1][y] && label2[x-1][y]>0)label2[x][y]=label2[x-1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y+1] && label2[x][y+1]>0)label2[x][y]=label2[x][y+1] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y-1] && label2[x][y-1]>0)label2[x][y]=label2[x][y-1] ;
	}

	for (i=120;i=0;i--) {
	int x = i / WIDTH;
	int y = i - (i / WIDTH)*WIDTH;    //
	if (label2[x][y]!=0 && label2[x][y]>label2[x+1][y] && label2[x+1][y]>0)label2[x][y]=label2[x+1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x-1][y] && label2[x-1][y]>0)label2[x][y]=label2[x-1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y+1] && label2[x][y+1]>0)label2[x][y]=label2[x][y+1] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y-1] && label2[x][y-1]>0)label2[x][y]=label2[x][y-1] ;
	}
	*/




	for (i = 0; i<BOARD_MAX; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // jinti

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x + 1][y]; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x - 1][y]; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y + 1]; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y - 1]; }

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}


	for (i = BOARD_MAX-1; i <1; i--) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    //

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x + 1][y]; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x - 1][y]; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y + 1]; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y - 1]; }

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}

	for (i = 0; i<BOARD_MAX; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // jinti

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x + 1][y]; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x - 1][y]; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y + 1]; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y - 1]; }

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}




	for (i = 0; i<BOARD_MAX; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // jinti

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}

		for (i = BOARD_MAX-1; i < 1; i--) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    //

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}


		
	int turn1 = 0;
	int kanou = 0;

	for (i = 0; i<BOARD_MAX; i++) {
		label_num[i] = 0;

	}

	for (i = 0; i<BOARD_MAX; i++) {
		if (label1[i - (i / WIDTH)*WIDTH][i / WIDTH]>0 && label1[i - (i / WIDTH)*WIDTH][i / WIDTH]<3) { turn1 = turn1 + 1; }
		if (label1[i - (i / WIDTH)*WIDTH][i / WIDTH] == 0 && jinti1[label2[i - (i / WIDTH)*WIDTH][i / WIDTH]]>3) { kanou = kanou + 1; };

		if (label1[i - (i / WIDTH)*WIDTH][i / WIDTH] == 0 && label2[i - (i / WIDTH)*WIDTH][i / WIDTH]>0)label_num[label2[i - (i / WIDTH)*WIDTH][i / WIDTH]] = label_num[label2[i - (i / WIDTH)*WIDTH][i / WIDTH]] + 1;

	}


  	prt("label=%d \n", label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]);
	prt("jinti=%d \n", jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]);
	prt("tz=%d \n", tz);

	for (i = 0; i<BOARD_MAX; i++) {
		if (label_num[i]>0)prt("label_num=%d %d \n", i, label_num[i]);
	}


	// label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]>0
	//  turn1>6 &&


	if (jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]>3)  return 0;
	//	if (jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]>7)  return 0;
	//	if (jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]<3)  return 6;

	if (turn1>15 && jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]<4)  return 5;
	if (label1[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]!=0)  return 6;

//	if ((tz - (tz / WIDTH)*WIDTH)==15 && (tz / WIDTH)==10)  return 6;


	//	if (turn1<11 && jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]<3)  return 5;
	//if (turn1>10 && jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]<7)  return 5;


	/*
	if (kanou<1) {  // pass
		return 7;
	}
		*/
//	return 7;
	return 0;

}


// put stone. success returns 0. in playout, fill_eye_err = 1.
int put_stone2(int tz, int color, int fill_eye_err)
{
	int around[4][3];
	int un_col = flip_color(color);
	int space = 0;
	int wall = 0;
	int mycol_safe = 0;
	int capture_sum = 0;
	int ko_maybe = 0;
	int liberty, stone;
	int i;

	int label1[13][13];
	int label2[13][13];
	int number = 1;
	int jinti1[121];

	if (tz == 0) {  // pass
		if (ko_z != 0) hash_xor(ko_z, HASH_KO);
		ko_z = 0;
		hash_pass();
		return 0;
	}

	// count 4 neighbor's liberty and stones.
	for (i = 0; i<4; i++) {
		int z, c, liberty, stone;
		around[i][0] = around[i][1] = around[i][2] = 0;
		z = tz + dir4[i];
		c = board[z];  // color
		if (c == 0) space++;
		if (c == 3) wall++;
		if (c == 0 || c == 3) continue;
		count_liberty(z, &liberty, &stone);
		around[i][0] = liberty;
		around[i][1] = stone;
		around[i][2] = c;
		if (c == un_col && liberty == 1) { capture_sum += stone; ko_maybe = z; }
		if (c == color  && liberty >= 2) mycol_safe++;
	}

	if (capture_sum == 0 && space == 0 && mycol_safe == 0) return 1; // suicide
	if (tz == ko_z) return 2; // ko
	if (wall + mycol_safe == 4 && fill_eye_err) return 3; // eye
	if (board[tz] != 0) return 4;



	//argon code 1 \83A\83\8B\83S\83\93\81@\83G\83N\83X\83g\83\89\81@\83R\81[\83h
	for (i = 0; i<121; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // 106 = 9*11 + 7 = (x,y)=(7,9) -> 79

		label2[x][y] = 0;
	}



	for (i = 0; i<121; i++) {
		jinti1[i] = 8;

	}


	for (i = 0; i<121; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // 106 = 9*11 + 7 = (x,y)=(7,9) -> 79
		label1[x][y] = board[i];
		//			        prt("%d ",label1[x][y]);

	}

	for (i = 0; i<121; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // labeling



		if (label1[x][y] == 0 && label2[x][y] == 0 && label2[x + 1][y] == 0 && label2[x - 1][y] == 0 && label2[x][y + 1] == 0 && label2[x][y - 1] == 0
			) {
			label2[x][y] = number;
			number = number + 1;
		}
		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] == label2[x - 1][y];
		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] == label2[x][y - 1];
		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] == label2[x + 1][y];
		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] == label2[x][y + 1];

	}


	for (i = 120; i = 0; i--) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // 

		if (label1[x][y] == 0 && label2[x][y]>label2[x - 1][y] && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] = label2[x - 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y - 1] && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] = label2[x][y - 1];
		if (label1[x][y] == 0 && label2[x][y]>label2[x + 1][y] && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] = label2[x + 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y + 1] && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] = label2[x][y + 1];
	}


	for (i = 0; i<121; i++) {
		int x = i / WIDTH;
		int y = i - (i / WIDTH)*WIDTH;    // 

		if (label1[x][y] == 0 && label2[x][y]>label2[x - 1][y] && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] = label2[x - 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y - 1] && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] = label2[x][y - 1];
		if (label1[x][y] == 0 && label2[x][y]>label2[x + 1][y] && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] = label2[x + 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y + 1] && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] = label2[x][y + 1];
	}


	for (i = 120; i = 0; i--) {
		int x = i / WIDTH;
		int y = i - (i / WIDTH)*WIDTH;    // 

		if (label1[x][y] == 0 && label2[x][y]>label2[x - 1][y] && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] = label2[x - 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y - 1] && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] = label2[x][y - 1];
		if (label1[x][y] == 0 && label2[x][y]>label2[x + 1][y] && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] = label2[x + 1][y];
		if (label1[x][y] == 0 && label2[x][y]>label2[x][y + 1] && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] = label2[x][y + 1];
	}

	for (i = 0; i<121; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // labeling

		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x - 1][y] == 0 && label2[x - 1][y]>0)label2[x][y] == label2[x - 1][y];
		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x][y - 1] == 0 && label2[x][y - 1]>0)label2[x][y] == label2[x][y - 1];
		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x + 1][y] == 0 && label2[x + 1][y]>0)label2[x][y] == label2[x + 1][y];
		if (label1[x][y] == 0 && label2[x][y] == 0 && label1[x][y + 1] == 0 && label2[x][y + 1]>0)label2[x][y] == label2[x][y + 1];

	}
	/*

	*/




	/*
	for (i=0;i<121;i++) {
	int x = i / WIDTH;
	int y = i - (i / WIDTH)*WIDTH;    //
	if (label2[x][y]!=0 && label2[x][y]>label2[x+1][y] && label2[x+1][y]>0)label2[x][y]=label2[x+1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x-1][y] && label2[x-1][y]>0)label2[x][y]=label2[x-1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y+1] && label2[x][y+1]>0)label2[x][y]=label2[x][y+1] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y-1] && label2[x][y-1]>0)label2[x][y]=label2[x][y-1] ;
	}

	for (i=120;i=0;i--) {
	int x = i / WIDTH;
	int y = i - (i / WIDTH)*WIDTH;    //
	if (label2[x][y]!=0 && label2[x][y]>label2[x+1][y] && label2[x+1][y]>0)label2[x][y]=label2[x+1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x-1][y] && label2[x-1][y]>0)label2[x][y]=label2[x-1][y] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y+1] && label2[x][y+1]>0)label2[x][y]=label2[x][y+1] ;
	if (label2[x][y]!=0 && label2[x][y]>label2[x][y-1] && label2[x][y-1]>0)label2[x][y]=label2[x][y-1] ;
	}
	*/



	for (i = 0; i<121; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // jinti

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x + 1][y]; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x - 1][y]; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y + 1]; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y - 1]; }

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}

	for (i = 120; i = 0; i--) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // 

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x + 1][y]; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x - 1][y]; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y + 1]; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y - 1]; }

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}

	for (i = 0; i<121; i++) {
		int y = i / WIDTH;
		int x = i - (i / WIDTH)*WIDTH;    // jinti

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x + 1][y]; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x - 1][y]; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y + 1]; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] == 8) { jinti1[label2[x][y]] = label1[x][y - 1]; }

		if (label1[x][y] == 0 && label1[x + 1][y]<3 && label1[x + 1][y]>0 && jinti1[label2[x][y]] != label1[x + 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x - 1][y]<3 && label1[x - 1][y]>0 && jinti1[label2[x][y]] != label1[x - 1][y]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y + 1]<3 && label1[x][y + 1]>0 && jinti1[label2[x][y]] != label1[x][y + 1]) { jinti1[label2[x][y]] = 9; }
		if (label1[x][y] == 0 && label1[x][y - 1]<3 && label1[x][y - 1]>0 && jinti1[label2[x][y]] != label1[x][y - 1]) { jinti1[label2[x][y]] = 9; }

	}


	int turn1 = 0;
	int kanou = 0;
	for (i = 0; i<121; i++) {
		if (label1[i - (i / WIDTH)*WIDTH][i / WIDTH]>0 && label1[i - (i / WIDTH)*WIDTH][i / WIDTH]<3) { turn1 = turn1 + 1; }
		if (label1[i - (i / WIDTH)*WIDTH][i / WIDTH] == 0 && jinti1[label2[i - (i / WIDTH)*WIDTH][i / WIDTH]]>7) { kanou = kanou + 1; };
	}

	//  prt("%d ",label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]);
	//      prt("%d ",jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]);

	if (turn1<11 && jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]<3)  return 5;
	if (turn1>10 && jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]<9)  return 5;

	// label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]>0
	//  turn1>6 &&

	if (turn1>10 && kanou<1) {  // pass
		return 7;
	}

	for (i = 0; i<4; i++) {
		int lib = around[i][0];
		int c = around[i][2];
		if (c == un_col && lib == 1 && board[tz + dir4[i]] != 0) {
			take_stone(tz + dir4[i], un_col);
		}
	}

	board[tz] = color;
	hash_xor(tz, color);
	hash_pass();
	if (ko_z != 0) hash_xor(ko_z, HASH_KO);

	count_liberty(tz, &liberty, &stone);
	if (capture_sum == 1 && stone == 1 && liberty == 1) {
		ko_z = ko_maybe;
		hash_xor(ko_z, HASH_KO);
	}
	else {
		ko_z = 0;
	}


	/*
	//pass
	int kanou=0;
	for (i=0;i<121;i++) {
	if (jinti1[label2[i - (i / WIDTH)*WIDTH][i / WIDTH]]>8){kanou=kanou+1;}
	}
	if (turn1>6 && kanou<50 ) {  // pass
	return 0;

	}
	*/
	//	if (turn1<11 && jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]>3)  return 0;
	if (jinti1[label2[tz - (tz / WIDTH)*WIDTH][tz / WIDTH]]>3)  return 0;
	return 0;
}



int count_score(int turn_color)
{
	int x, y, i;
	int score = 0, win;
	int black_area = 0, white_area = 0, black_sum, white_sum;
	int mk[4];
	int kind[3];

	kind[0] = kind[1] = kind[2] = 0;
	for (y = 0; y<B_SIZE; y++) for (x = 0; x<B_SIZE; x++) {
		int z = get_z(x + 1, y + 1);
		int c = board[z];
		kind[c]++;
		if (c != 0) {
			if (c == 1) board_area_sum[z]++;
			if (c == 2) board_area_sum[z]--;
			continue;
		}
		mk[1] = mk[2] = 0;
		for (i = 0; i<4; i++) mk[board[z + dir4[i]]]++;
		if (mk[1] && mk[2] == 0) {
			black_area++;
			board_area_sum[z]++;
		}
		if (mk[2] && mk[1] == 0) {
			white_area++;
			board_area_sum[z]--;
		}
	}

	black_sum = kind[1] + black_area;
	white_sum = kind[2] + white_area;
	score = black_sum - white_sum;

	win = 0;
	if (score - komi > 0) win = 1;


	if (win == 1) { // black win
		for (y = 0; y<B_SIZE; y++) for (x = 0; x<B_SIZE; x++) {
			int z = get_z(x + 1, y + 1);
			if (board[z] == 1) board_winner[0][z]++;
		}
		winner_count[0]++;
	}
	else {
		for (y = 0; y<B_SIZE; y++) for (x = 0; x<B_SIZE; x++) {
			int z = get_z(x + 1, y + 1);
			if (board[z] == 2) board_winner[1][z]++;
		}
		winner_count[1]++;
	}



	if (turn_color == 2) win = -win;

	//prt("black_sum=%2d, (stones=%2d, area=%2d)\n",black_sum, kind[1], black_area);
	//prt("white_sum=%2d, (stones=%2d, area=%2d)\n",white_sum, kind[2], white_area);
	//prt("score=%d, win=%d\n",score, win);
	return win;
}




char *pattern3x3[] = {
	"X.X"
	"..."
	"XXX",

	".XX"
	"X.."
	"XO.",

	".XO"
	"X.."
	"O..",

	"X.."
	"..X"
	"...",

	".X."
	"X.O"
	".O.",

	"OXO"
	"X.."
	"O..",


	".X."
	"X.X"
	"..X",

	"X.."
	"..X"
	"...",


	".XO"
	"X.."
	"O.X",

	"..."
	"..X"
	".OX",

	"OX."
	"X.O"
	"...",

	"X.."
	"..X"
	"...",

	"XXO"
	"X.X"
	"OX.",

	"OX."
	"X.."
	"...",

	"..."
	"..X"
	"..O",

	".OX"
	"O.O"
	"XXO",

	"..."
	"..X"
	".XO",

	"..O"
	"..."
	"OX.",


	NULL
};

#define EXPAND_PATTERN_MAX  (8*500)
int e_pat_num = 0;
int e_pat[EXPAND_PATTERN_MAX][9];      // rotate and flip pattern
int e_pat_bit[EXPAND_PATTERN_MAX][2];  // [0] ...pattern, [1]...mask
int dir_3x3[9] = { -WIDTH - 1, -WIDTH, -WIDTH + 1, -1, 0, +1, +WIDTH - 1, +WIDTH, +WIDTH + 1 };

void expand_pattern3x3()
{
	int n, i, j;
	e_pat_num = 0;
	for (n = 0;; n++) {
		int i, j;
		char *p = pattern3x3[n];
		if (p == NULL) break;
		if (e_pat_num > EXPAND_PATTERN_MAX - 8) { prt("e_pat_num over Err\n"); exit(0); }

		for (i = 0; i<9; i++) {
			int m = 0;
			char c = *(p + i);
			if (c == '.') m = 0;
			if (c == 'X') m = 1;
			if (c == 'O') m = 2;
			if (c == '#') m = 3;
			if (c == '?') m = 4;
			e_pat[e_pat_num][i] = m;
		}
		e_pat_num++;
		for (i = 0; i<2; i++) {
			int *p;
			int *q;
			for (j = 0; j<3; j++) {
				p = e_pat[e_pat_num - 1];
				q = e_pat[e_pat_num];
				// roteta 90
				//  "012"      "630"
				//  "345"  ->  "741"
				//  "678"      "852"
				q[0] = p[6];
				q[1] = p[3];
				q[2] = p[0];
				q[3] = p[7];
				q[4] = p[4];
				q[5] = p[1];
				q[6] = p[8];
				q[7] = p[5];
				q[8] = p[2];
				e_pat_num++;
			}
			if (i == 1) break;
			p = e_pat[e_pat_num - 1];
			q = e_pat[e_pat_num];
			// vertical flip
			q[0] = p[6];
			q[1] = p[7];
			q[2] = p[8];
			q[3] = p[3];
			q[4] = p[4];
			q[5] = p[5];
			q[6] = p[0];
			q[7] = p[1];
			q[8] = p[2];
			e_pat_num++;
		}
	}

	for (i = 0; i<e_pat_num; i++) {
		e_pat_bit[i][0] = 0;
		e_pat_bit[i][1] = 0;
		//  prt("%4d\n",i);
		for (j = 0; j<9; j++) {
			int c = e_pat[i][j];
			int mask = 3;
			//    prt("%d",c);
			//    if ((j+1)%3==0) prt("\n");
			if (c == 4) {
				mask = 0;
				c = 0;
			}
			e_pat_bit[i][0] = e_pat_bit[i][0] << 2;
			e_pat_bit[i][1] = e_pat_bit[i][1] << 2;
			e_pat_bit[i][0] |= c;
			e_pat_bit[i][1] |= mask;
		}
		//  prt("bit=%08x,mask=%08x\n",e_pat_bit[i][0],e_pat_bit[i][1]);
	}
	prt("pattern3x3 num=%d, e_pat_num=%d\n", n, e_pat_num);
}

// return pattern number, -1 ... not found 
int match_pattern3x3(int z, int col)
{
#if 1    // 413 playouts/sec
	int pat_bit = 0;
	int i, j;
	for (j = 0; j<9; j++) {
		int c = board[z + dir_3x3[j]];
		if (col == 2 && (c == 1 || c == 2)) c = 3 - c;
		pat_bit = pat_bit << 2;
		pat_bit |= c;
	}
	for (i = 0; i<e_pat_num; i++) {
		int e_bit = e_pat_bit[i][0];
		int e_mask = e_pat_bit[i][1];
		if (e_bit == (pat_bit & e_mask)) return i;
	}
	return -1;
#else    // 353 playouts/sec
	int i, j;
	for (i = 0; i<e_pat_num; i++) {
		for (j = 0; j<9; j++) {
			int c = board[z + dir_3x3[j]];
			int e = e_pat[i][j];
			if (col == 2 && e == 1) e = 2;
			else if (col == 2 && e == 2) e = 1;
			if (e == 4) continue;
			if (c != e) break;
		}
		if (j == 9) return i;
	}
	return -1;
#endif
}


int get_prob(int z, int prev_z, int col)
{
	const int MAX_PROB = INT_MAX / BOARD_MAX;
	int pr = 100;
	int un_col = flip_color(col);
	int y = z / WIDTH;
	int x = z - y*WIDTH;
	int prev_y = prev_z / WIDTH;
	int prev_x = prev_z - prev_y*WIDTH;
	int dx = abs(x - prev_x);
	int dy = abs(y - prev_y);
	int m;
	int i;
	int sc[4];

	// add your code, start -->


	sc[0] = sc[1] = sc[2] = sc[3] = 0;
	for (i = 0; i<8; i++) {
		sc[board[z + dir8[i]]]++;
	}
	if (sc[un_col] >= 3 && sc[col] == 0) pr /= 2;
	if (sc[col] >= 3 && sc[un_col] == 0) pr *= 2;

	m = -1;
	if (prev_z != 0 && ((dx + dy) == 1 || (dx*dy) == 1)) {
		m = match_pattern3x3(z, col);
	}
	if (m >= 0) {
		int n = m / 18;  // pattern number
						 //  prt("match=%3d,z=%s,col=%d\n",m,get_char_z(z),col);
		pr *= 1000;
		if (n == 0) pr *= 2;
	}


	// add your code, end <--

	if (pr < 1) pr = 1;
	if (pr > MAX_PROB) pr = MAX_PROB;
	return pr;
}





int playout(int turn_color)
{
	int color = turn_color;
	int previous_z = 0;
	int loop;
	int loop_max = B_SIZE*B_SIZE + 900;  // for triple ko


/*
		//use caffe model 
			int x3,y3;

			prt("\n");
			prt("polnet1");
			prt("\n");

		for (y3 = 0; y3<B_SIZE; y3++) {
			for (x3 = 0; x3<B_SIZE; x3++) {
//			int z = get_z(x + 1, y + 1);
			prt("%5.3f ", board3[xy2pos(x3,y3)]);
		}
			prt("\n");

	}
*/


	all_playouts++;
	for (loop = 0; loop<loop_max; loop++) {
		// all empty points are candidates.
		int empty[BOARD_MAX][2];  // [0]...z, [1]...probability
		int empty_num = 0;
		int prob_sum = 0;
		int x, y, z, err, pr,b1,c1;
		for (y = 0; y<B_SIZE; y++) for (x = 0; x<B_SIZE; x++) {
			int z = get_z(x + 1, y + 1);
			if (board[z] != 0) continue;
			empty[empty_num][0] = z;

			//pr = 1;
			//pr = get_prob(z, previous_z, color);
//			int polnet(int color);
//			polnet(color);


			//testCNN(color);
			//int pr;
			float a1;
			b1=1;
			a1=board3[xy2pos(x,y)];
			//prt("%5.3f ", a1);
			b1=a1;
			c1=b1;
			if (b1>0) pr=c1;
			if (b1<1)pr=1;

			empty[empty_num][1] = pr;
			prob_sum += pr;
			empty_num++;
		}
		for (;;) {
			int i = 0;
			if (empty_num == 0) {
				z = 0;
			}
			else {
				int r = rand() % prob_sum;
				int sum = 0;
				for (i = 0; i<empty_num; i++) {
					sum += empty[i][1];    // 0,1,2   [0]=1, [1]=1, [2]=1 
					if (sum > r) break;
				}
				if (i == empty_num) { prt("Err! prob_sum=%d,sum=%d,r=%d,r=%d\n", prob_sum, sum, r, i); exit(0); }
				z = empty[i][0];
			}
			err = put_stone(z, color, FILL_EYE_ERR);
			if (err == 0) break;  // pass is ok.
			prob_sum -= empty[i][1];
			empty[i][0] = empty[empty_num - 1][0];  // err, delete
			empty[i][1] = empty[empty_num - 1][1];
			empty_num--;
		}
		if (flag_test_playout) record[moves++] = z;

		if (depth < D_MAX) path2[depth++] = z;

		if (z == 0 && previous_z == 0) break;  // continuous pass
		previous_z = z;
		//  prt("loop=%d,z=%s,c=%d,empty_num=%d,ko_z=%d\n",loop,get_char_z(z),color,empty_num,ko_z);
		color = flip_color(color);
	}
	return count_score(turn_color);
}

int primitive_monte_calro(int color)
{
	int    try_num = 30; // number of playout
	int    best_z = 0;
	double best_value;
	double win_rate;
	int x, y, err, i, win_sum, win;

	int ko_z_copy;
	int board_copy[BOARD_MAX];  // keep current board
	ko_z_copy = ko_z;
	memcpy(board_copy, board, sizeof(board));

	best_value = -100;

	// try all empty point
	for (y = 0; y<B_SIZE; y++) for (x = 0; x<B_SIZE; x++) {
		int z = get_z(x + 1, y + 1);
		if (board[z] != 0) continue;

		err = put_stone(z, color, FILL_EYE_ERR);
		if (err != 0) continue;

		win_sum = 0;
		for (i = 0; i<try_num; i++) {
			int board_copy2[BOARD_MAX];
			int ko_z_copy2 = ko_z;
			memcpy(board_copy2, board, sizeof(board));

			win = -playout(flip_color(color));
			win_sum += win;

			ko_z = ko_z_copy2;
			memcpy(board, board_copy2, sizeof(board));
		}
		win_rate = (double)win_sum / try_num;
		//  print_board();
		//  prt("z=%d,win=%5.3f\n",get81(z),win_rate);

		if (win_rate > best_value) {
			best_value = win_rate;

			best_z = z;

			//    prt("best_z=%d,color=%d,v=%5.3f,try_num=%d\n",get81(best_z),color,best_value,try_num);
		}

		ko_z = ko_z_copy;
		memcpy(board, board_copy, sizeof(board));  // resume board
	}
	return best_z;
}



// following are for UCT

typedef struct {
	int    z;          // move position
	int    games;      // number of games
	double rate;       // winrate
	int    rave_games; // (RAVE) number of games
	double rave_rate;  // (RAVE) winrate
	int    next;       // next node
	double bonus;      // shape bonus
} CHILD;

#define CHILD_SIZE  (B_SIZE*B_SIZE+1)  // +1 for PASS

typedef struct {
	int child_num;
	CHILD child[CHILD_SIZE];
	int child_games_sum;
	int child_rave_games_sum;
} NODE;

#define NODE_MAX 90000
NODE node[NODE_MAX];
int node_num = 0;
const int NODE_EMPTY = -1; // no next node
const int ILLEGAL_Z = -1; // illegal move


void add_child(NODE *pN, int z, double bonus)
{
	int n = pN->child_num;
	pN->child[n].z = z;
	pN->child[n].games = 0;
	pN->child[n].rate = 0;
	pN->child[n].rave_games = 0;
	pN->child[n].rave_rate = 0;
	pN->child[n].next = NODE_EMPTY;
	pN->child[n].bonus = bonus;  // from 0 to 10, good move has big bonus.
	pN->child_num++;
}

// create new node. return node index.
int create_node(int prev_z)
{
	int x, y, z, i, j;
	NODE *pN;

	if (node_num == NODE_MAX) { prt("node over Err\n"); exit(0); }
	pN = &node[node_num];
	pN->child_num = 0;
	pN->child_games_sum = 0;
	pN->child_rave_games_sum = 0;

	for (y = 0; y<B_SIZE; y++) for (x = 0; x<B_SIZE; x++) {
		z = get_z(x + 1, y + 1);
		if (board[z] != 0) continue;
		add_child(pN, z, 0);
	}
	add_child(pN, 0, 0);  // add PASS

						  // sort children
	for (i = 0; i<pN->child_num - 1; i++) {
		double max_b = pN->child[i].bonus;
		int    max_i = i;
		CHILD tmp;
		for (j = i + 1; j<pN->child_num; j++) {
			CHILD *c = &pN->child[j];
			if (max_b >= c->bonus) continue;
			max_b = c->bonus;
			max_i = j;
		}
		if (max_i == i) continue;
		tmp = pN->child[i];
		pN->child[i] = pN->child[max_i];
		pN->child[max_i] = tmp;
	}

	node_num++;
	return node_num - 1;
}

int select_best_ucb(int node_n, int color)
{
	NODE *pN = &node[node_n];
	int select = -1;
	double max_ucb = -999;
	double ucb = 0, ucb_rave = 0, beta;
	int i;

	for (i = 0; i<pN->child_num; i++) {
		CHILD *c = &pN->child[i];
		if (c->z == ILLEGAL_Z) continue;

		if (c->games == 0) {
			ucb_rave = 10000 + (rand() & 0x7fff);  // try once
		}
		else {
			const double C = 0.30;    // depends on program
			const double RAVE_D = 3000;
			double moveCount = c->games;
			double raveCount = c->rave_games;
			double rave = c->rave_rate;
			if (c->z == 0) {  // dont select pass
				rave = 1 - color;
				raveCount = pN->child_games_sum;
			}

			beta = raveCount / (raveCount + moveCount + raveCount*moveCount / RAVE_D);

			ucb = c->rate + C * sqrt(log((double)pN->child_games_sum) / c->games);

			ucb_rave = beta * rave + (1 - beta) * ucb;
			//    if ( depth==0 ) prt("%2d:z=%2d,rate=%6.3f,games=%4d, rave_r=%6.3f,g=%4d, beta=%f,ucb_rave=%f\n", i, get81(c->z), c->rate, c->games, c->rave_rate, c->rave_games,beta,ucb_rave);
		}
		if (ucb_rave > max_ucb) {
			max_ucb = ucb_rave;
			select = i;
		}
	}
	if (select == -1) { prt("Err! select\n"); exit(0); }
	return select;
}

void update_rave(NODE *pN, int color, int current_depth, double win)
{
	int played_color[BOARD_MAX];
	int i, z;
	int c = color;

	memset(played_color, 0, sizeof(played_color));
	for (i = current_depth; i<depth; i++) {
		z = path2[i];
		if (played_color[z] == 0) played_color[z] = c;
		c = flip_color(c);
	}

	played_color[0] = 0;	// ignore pass

	for (i = 0; i<pN->child_num; i++) {
		CHILD *c = &pN->child[i];
		if (c->z == ILLEGAL_Z) continue;
		if (played_color[c->z] != color) continue;
		c->rave_rate = (c->rave_games * c->rave_rate + win) / (c->rave_games + 1);
		c->rave_games++;
		pN->child_rave_games_sum++;
	}
}

int search_uct(int color, int node_n)
{
	NODE *pN = &node[node_n];
	CHILD *c = NULL;
	int select, z, err, win, current_depth;
	for (;;) {
		select = select_best_ucb(node_n, color);
		c = &pN->child[select];
		z = c->z;
		err = put_stone(z, color, FILL_EYE_ERR);
		if (err == 0) break;

		/*
		//argon extra code
		if (err == 7) {
			z = 0;
			break;
		}
		*/
		c->z = ILLEGAL_Z;     // select other move
	}

	current_depth = depth;
	path2[depth++] = c->z;

	// playout in first time. <= 10 can reduce node.
	if (c->games <= 0 || depth == D_MAX || (c->z == 0 && depth >= 2 && path2[depth - 2] == 0)) {
		win = -playout(flip_color(color));
	}
	else {
		if (c->next == NODE_EMPTY) c->next = create_node(c->z);
		win = -search_uct(flip_color(color), c->next);
	}

	update_rave(pN, color, current_depth, win);

	// update winrate
	c->rate = (c->rate * c->games + win) / (c->games + 1);
	c->games++;
	pN->child_games_sum++;
	return win;
}


// get mill second time. clock() returns process CPU times on Linux, not proper when multi thread.
double get_clock()
{
#if defined(_MSC_VER)
	return clock();
#else
	struct timeval  val;
	struct timezone zone;
	if (gettimeofday(&val, &zone) == -1) { prt("time err\n"); exit(0); }
	double t = val.tv_sec*1000.0 + (val.tv_usec / 1000.0);
	return t;
#endif
}

// get sec time.
double get_spend_time(double ct)
{
	//int div = CLOCKS_PER_SEC;	// 1000 ...VC, 1000000 ...gcc
	int div = 1000;
	return (double)(get_clock() + 1 - ct) / div;
}

double start_time;
double time_limit_sec = 3.0;

int is_time_over()
{
	if (use_time_control == 0) return 0;
	if (get_spend_time(start_time) >= time_limit_sec) return 1;
	return 0;
}

double count_total_time()
{
	int i;
	double total_time[2];

	total_time[0] = 0;  // black time
	total_time[1] = 0;  // white

	for (i = 0; i<moves; i++) {
		total_time[i & 1] += record_time[i];
	}
	return total_time[moves & 1];
}


int uct_loop = 65000;  // number of uct loop

int get_best_uct(int color)
{
	int next, i, best_z, best_i = -1;
	int max = -999;
	NODE *pN;
	int prev_z = 0;

	int first = 0;
	int second = 0;
	int third = 0;
	int recent = 0;
	int temp4 = 0;
	int temp5 = 0;
	int temp6 = 0;
	int temp7 = 0;
	int temp8 = 0;
	int temp9 = 0;
	int temp10 = 0;
	int temp11 = 0;
	int temp12 = 0;
	int temp13 = 0;
	int temp14 = 0;
	int temp15 = 0;
	int temp16 = 0;
	int temp17 = 0;
	int temp18 = 0;
	int temp19 = 0;
	int temp20 = 0;
	int temp21 = 0;
	int temp22 = 0;
	int temp23 = 0;
	int temp24 = 0;
	int temp25 = 0;



	int tempi1 = 0;
	int tempi2 = 0;
	int tempi3 = 0;
	int tempi4 = 0;
	int tempi5 = 0;
	int tempi6 = 0;
	int tempi7 = 0;
	int tempi8 = 0;
	int tempi9 = 0;
	int tempi10 = 0;
	int tempi11 = 0;
	int tempi12 = 0;
	int tempi13 = 0;
	int tempi14 = 0;
	int tempi15 = 0;
	int tempi16 = 0;
	int tempi17 = 0;
	int tempi18 = 0;
	int tempi19 = 0;
	int tempi20 = 0;
	int tempi21 = 0;
	int tempi22 = 0;
	int tempi23 = 0;
	int tempi24 = 0;
	int tempi25 = 0;



	int temp1_z;
	int temp2_z;
	int temp3_z;
	int temp4_z;
	int temp5_z;
	int temp6_z;
	int temp7_z;
	int temp8_z;
	int temp9_z;
	int temp10_z;
	int temp11_z;
	int temp12_z;
	int temp13_z;
	int temp14_z;
	int temp15_z;
	int temp16_z;
	int temp17_z;
	int temp18_z;
	int temp19_z;
	int temp20_z;
	int temp21_z;
	int temp22_z;
	int temp23_z;
	int temp24_z;
	int temp25_z;


	if (moves > 0) prev_z = record[moves - 1];
	node_num = 0;
	next = create_node(prev_z);




		//use caffe model 
			int x3,y3;

		for (y3 = 0; y3<B_SIZE; y3++) {
			for (x3 = 0; x3<B_SIZE; x3++) {
			int z = get_z(x3 + 1, y3 + 1);
			board2[xy2pos(x3,y3)]=board[z];
		}
	}


			//polycy net
			polnet(color);

			//testCNN(color);

			prt("\n");
			prt("polnet0");
			prt("\n");

		for (y3 = 0; y3<B_SIZE; y3++) {
			for (x3 = 0; x3<B_SIZE; x3++) {
//			int z = get_z(x + 1, y + 1);
			prt("%5.3f ", board3[xy2pos(x3,y3)]);
		}
			prt("\n");

	}


	for (i = 0; i<uct_loop; i++) {
		int board_copy[BOARD_MAX];
		int ko_z_copy = ko_z;
		memcpy(board_copy, board, sizeof(board));

		depth = 0;
		search_uct(color, next);

		ko_z = ko_z_copy;
		memcpy(board, board_copy, sizeof(board));

		if (is_time_over()) break;
	}
	pN = &node[next];
	for (i = 0; i<pN->child_num; i++) {
		CHILD *c = &pN->child[i];


		if (c->games > max) {


			tempi25 = tempi24;
			tempi24 = tempi23;
			tempi23 = tempi22;
			tempi22 = tempi21;
			tempi21 = tempi20;
			tempi20 = tempi19;
			tempi19 = tempi18;
			tempi18 = tempi17;
			tempi17 = tempi16;
			tempi16 = tempi15;
			tempi15 = tempi14;
			tempi14 = tempi13;
			tempi13 = tempi12;
			tempi12 = tempi11;
			tempi11 = tempi10;
			tempi10 = tempi9;
			tempi9 = tempi8;
			tempi8 = tempi7;
			tempi7 = tempi6;
			tempi6 = tempi5;
			tempi5 = tempi4;
			tempi4 = tempi3;
			tempi3 = tempi2;
			tempi2 = tempi1;

			temp25 = temp24;
			temp24 = temp23;
			temp23 = temp22;
			temp22 = temp21;
			temp21 = temp20;
			temp20 = temp19;
			temp19 = temp18;
			temp18 = temp17;
			temp17 = temp16;
			temp16 = temp15;
			temp15 = temp14;
			temp14 = temp13;
			temp13 = temp12;
			temp12 = temp11;
			temp11 = temp10;
			temp10 = temp9;
			temp9 = temp8;
			temp8 = temp7;
			temp7 = temp6;
			temp6 = temp5;
			temp5 = temp4;
			temp4 = third;
			third = second;
			second = first;

			best_i = i;
			tempi1 = i;
			max = c->games;
			first = c->games;

		}
		else if (c->games > 5) {

			temp25 = temp24;
			temp24 = temp23;
			temp23 = temp22;
			temp22 = temp21;
			temp21 = temp20;
			temp20 = temp19;
			temp19 = temp18;
			temp18 = temp17;
			temp17 = temp16;
			temp16 = temp15;

			temp15 = temp14;
			temp14 = temp13;
			temp13 = temp12;
			temp12 = temp11;
			temp11 = temp10;
			temp10 = temp9;
			temp9 = temp8;
			temp8 = temp7;
			temp7 = temp6;
			temp6 = temp5;
			temp5 = temp4;
			temp4 = third;
			third = second;
			second = first;


			tempi25 = tempi24;
			tempi24 = tempi23;
			tempi23 = tempi22;
			tempi22 = tempi21;
			tempi21 = tempi20;
			tempi20 = tempi19;
			tempi19 = tempi18;
			tempi18 = tempi17;
			tempi17 = tempi16;
			tempi16 = tempi15;

			tempi15 = tempi14;
			tempi14 = tempi13;
			tempi13 = tempi12;
			tempi12 = tempi11;
			tempi11 = tempi10;
			tempi10 = tempi9;
			tempi9 = tempi8;
			tempi8 = tempi7;
			tempi7 = tempi6;
			tempi6 = tempi5;
			tempi5 = tempi4;
			tempi4 = tempi3;
			tempi3 = tempi2;
			tempi2 = tempi1;

			tempi1 = i;
			first = c->games;

		}


		prt("%2d:z=%2d,rate=%6.3f,games=%4d, rave_r=%6.3f,g=%4d\n",
			i, get81(c->z), c->rate, c->games, c->rave_rate, c->rave_games);
	}


	//#define NUM_DATA 15                           /* \83\\81[\83g\82\B7\82\E9\83f\81[\83^\82̐\94 */
	int x[] = { first, second, third, temp4, temp5, temp6, temp7, temp8, temp9, temp10,temp11,temp12,temp13,temp14,temp15,temp16,temp17,temp18,temp19,temp20,temp21,temp22,temp23,temp24,temp25 };
	int y[] = { tempi1,tempi2,tempi3, tempi4, tempi5, tempi6, tempi7, tempi8, tempi9,tempi10, tempi11,tempi12,tempi13,tempi14,tempi15,tempi16,tempi17,tempi18,tempi19,tempi20,tempi21,tempi22,tempi23,tempi24,tempi25 };




	/* \83o\83u\83\8B\83\\81[\83g\82\F0\8Ds\82\A4 */
	//	int BubSort(int x[],int y[], int n);

	{
		int i, j, temp, tempb, n;
		n = 25;
		for (i = 0; i < n - 1; i++) {
			for (j = n - 1; j > i; j--) {
				if (x[j - 1] < x[j]) {
					temp = x[j];
					tempb = y[j];
					x[j] = x[j - 1];
					y[j] = y[j - 1];
					x[j - 1] = temp;
					y[j - 1] = tempb;

				}
			}
		}
	}

	
	int goukei=0;
	int kakuritu=0;
	int tempc=0;
	int tempd=0;
	goukei=(x[0]+x[1]);

	srand((unsigned)time(NULL));
	kakuritu=(rand()%goukei+1);
	
	if (kakuritu<x[1]){
	tempc=x[0];
	x[0]=x[1];
	x[1]=tempc;

	tempd=y[0];
	y[0]=y[1];
	y[1]=tempd;
	}


	prt("testx=%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24]);
	prt("testy=%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11], y[12], y[13], y[14], y[15], y[16], y[17], y[18], y[19], y[20], y[21], y[22], y[23], y[24]);
//	prt("testy=%d,best)


	best_z = pN->child[best_i].z;
	temp1_z = pN->child[y[0]].z;
	temp2_z = pN->child[y[1]].z;
	temp3_z = pN->child[y[2]].z;
	temp4_z = pN->child[y[3]].z;
	temp5_z = pN->child[y[4]].z;
	temp6_z = pN->child[y[5]].z;
	temp7_z = pN->child[y[6]].z;
	temp8_z = pN->child[y[7]].z;
	temp9_z = pN->child[y[8]].z;
	temp10_z = pN->child[y[9]].z;
	temp11_z = pN->child[y[10]].z;
	temp12_z = pN->child[y[11]].z;
	temp13_z = pN->child[y[12]].z;
	temp14_z = pN->child[y[13]].z;
	temp15_z = pN->child[y[14]].z;
	temp16_z = pN->child[y[15]].z;
	temp17_z = pN->child[y[16]].z;
	temp18_z = pN->child[y[17]].z;
	temp19_z = pN->child[y[18]].z;
	temp20_z = pN->child[y[19]].z;
	temp21_z = pN->child[y[20]].z;
	temp22_z = pN->child[y[21]].z;
	temp23_z = pN->child[y[22]].z;
	temp24_z = pN->child[y[23]].z;
	temp25_z = pN->child[y[24]].z;


	prt("temp_z=%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", temp1_z, temp2_z, temp3_z, temp4_z, temp5_z, temp6_z, temp7_z, temp8_z, temp9_z,temp10_z, temp11_z, temp12_z, temp13_z, temp14_z, temp15_z);

	prt("best_z=%d,rate=%6.3f,games=%4d,playouts=%d,nodes=%d\n",
		get81(best_z), pN->child[best_i].rate, max, all_playouts, node_num);

//	return best_z;

	if (argon_check1(best_z) == 0) {
//		return best_z;
		return temp1_z;
	}

	else if (argon_check1(temp1_z) == 0 ) {
		return temp1_z;
	}
	else if (argon_check1(temp2_z) == 0 ) {
		return temp2_z;
	}
	else if (argon_check1(temp3_z) == 0 ) {
		return temp3_z;
	}
	else if (argon_check1(temp4_z) == 0 ) {
		return temp4_z;
	}
	else if (argon_check1(temp5_z) == 0 ) {
		return temp5_z;
	}
	else if (argon_check1(temp6_z) == 0 ) {
		return temp6_z;
	}
	else if (argon_check1(temp7_z) == 0 ) {
		return temp7_z;
	}
	else if (argon_check1(temp8_z) == 0 ) {
		return temp8_z;
	}
	else if (argon_check1(temp9_z) == 0 ) {
		return temp9_z;
	}
	else if (argon_check1(temp10_z) == 0) {
		return temp10_z;
	}
	else if (argon_check1(temp11_z) == 0) {
		return temp11_z;
	}
	else if (argon_check1(temp12_z) == 0) {
		return temp12_z;
	}
	else if (argon_check1(temp13_z) == 0) {
		return temp13_z;
	}
	else if (argon_check1(temp14_z) == 0) {
		return temp14_z;
	}
	else if (argon_check1(temp15_z) == 0) {
		return temp15_z;
	}
	else if (argon_check1(temp16_z) == 0) {
		return temp16_z;
	}
	else if (argon_check1(temp17_z) == 0) {
		return temp17_z;
	}
	else if (argon_check1(temp18_z) == 0) {
		return temp18_z;
	}
	else if (argon_check1(temp19_z) == 0) {
		return temp19_z;
	}
	else if (argon_check1(temp20_z) == 0) {
		return temp20_z;
	}
	else if (argon_check1(temp21_z) == 0) {
		return temp21_z;
	}
	else if (argon_check1(temp22_z) == 0) {
		return temp22_z;
	}
	else if (argon_check1(temp23_z) == 0) {
		return temp23_z;
	}
	else if (argon_check1(temp24_z) == 0) {
		return temp24_z;
	}
	else if (argon_check1(temp25_z) == 0) {
		return temp25_z;
	}

	else {
//		return best_z;
		return 0;
	}

}

void init_board()
{
	int i, x, y;
	for (i = 0; i<BOARD_MAX; i++) board[i] = 3;
	for (y = 0; y<B_SIZE; y++) for (x = 0; x<B_SIZE; x++) board[get_z(x + 1, y + 1)] = 0;
	moves = 0;
	ko_z = 0;
	hashcode = 0;
}

void add_moves(int z, int color, double sec)
{
	int err = put_stone(z, color, FILL_EYE_OK);
	if (err != 0) { prt("put_stone err=%d\n", err); exit(0); }
	record[moves] = z;
	record_time[moves] = sec;
	moves++;
	print_board();
	prt("hashcode="); prt_code64(hashcode); prt("\n");
}

const int SEARCH_PRIMITIVE = 0;
const int SEARCH_UCT = 1;

int play_computer_move(int color, int search)
{
	double sec;
	int z;

	double total_time = count_total_time();

	double base_time = 60 * 15;    // 10 minutes
	double left_time = base_time - total_time;
	int div = 70; // 12_9*9 40 ... 13x13, 70 ... 19x19

//	if (left_time<600)int div = 20;
//	if (left_time<300)int div = 12;

		time_limit_sec = left_time / div;

	//  if ( left_time > 500       ) time_limit_sec = 10.0;

//	if (left_time < 120) time_limit_sec = 10.0;
//	if (left_time < 60) time_limit_sec = 5.0;
	if (left_time < 30) time_limit_sec = 1.0;
	if (left_time < 10) time_limit_sec = 0.2;
	if (use_time_control != 0) {
		prt("time_limit_sec=%.1f, total=%.1f, left=%.1f\n", time_limit_sec, total_time, left_time);
	}
	start_time = get_clock();

	all_playouts = 0;
	memset(board_area_sum, 0, sizeof(board_area_sum));
	memset(board_winner, 0, sizeof(board_winner));
	memset(winner_count, 0, sizeof(winner_count));

	if (search == SEARCH_UCT) {
		z = get_best_uct(color);
	}
	else {
		z = primitive_monte_calro(color);
	}

	//print_board_area();
	//print_criticality();

	sec = get_spend_time(start_time);
	prt("z=%s,color=%d,moves=%d,playouts=%d, %.1f sec(%.0f po/sec),depth=%d\n",
		get_char_z(z), color, moves, all_playouts, sec, all_playouts / sec, depth);

	add_moves(z, color, sec);
	return z;
}

void undo()
{
	int moves_copy = moves - 1;
	int color = 1;
	int i;

	if (moves == 0) return;
	init_board();

	for (i = 0; i<moves_copy; i++) {
		int z = record[i];
		int err = put_stone(z, color, FILL_EYE_OK);
		if (err != 0) { prt("put_stone err=%d\n", err); exit(0); }
		color = flip_color(color);
	}
	moves = moves_copy;
	print_board();
}

// print SGF game record
void print_sgf()
{
	int i;
	prt("(;GM[1]SZ[%d]KM[%.1f]PB[]PW[]\n", B_SIZE, komi);
	for (i = 0; i<moves; i++) {
		int z = record[i];
		int y = z / WIDTH;
		int x = z - y*WIDTH;
		const char *sStone[2] = { "B", "W" };
		prt(";%s", sStone[i & 1]);
		if (z == 0) {
			prt("[]");
		}
		else {
			prt("[%c%c]", x + 'a' - 1, y + 'a' - 1);
		}
		if (((i + 1) % 10) == 0) prt("\n");
	}
	prt(")\n");
}

void selfplay()
{
	int color = 1;
	int z, search;

	for (;;) {
		if (color == 1) {
			search = SEARCH_UCT; //SEARCH_PRIMITIVE;
		}
		else {
			search = SEARCH_UCT;
		}
		z = play_computer_move(color, search);
		if (z == 0 && moves > 1 && record[moves - 2] == 0) break;
		if (moves > 600) break;  // too long
		color = flip_color(color);
	}

	print_sgf();
}

void test_playout()
{
	flag_test_playout = 1;
	playout(1);
	print_board();
	print_sgf();
}


#define STR_MAX 256
#define TOKEN_MAX 3

void gtp_loop()
{
  char str[STR_MAX];
  char sa[TOKEN_MAX][STR_MAX];
  char seps[] = " ";
  char *token;
  int x,y,z,ax, count;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);
  for (;;) {
    if ( fgets(str, STR_MAX, stdin)==NULL ) break;
//  prt("gtp<-%s",str);
    count = 0;
    token = strtok( str, seps );
    while ( token != NULL ) {
      strcpy(sa[count], token);
      count++;
      if ( count == TOKEN_MAX ) break;
      token = strtok( NULL, seps );
    }

    if ( strstr(sa[0],"boardsize")     ) {
    int new_board_size = atoi( sa[1] );
      send_gtp("= \n\n");
    } else if ( strstr(sa[0],"clear_board")   ) {
      init_board();
      send_gtp("= \n\n");
    } else if ( strstr(sa[0],"quit") ) {
      break;
    } else if ( strstr(sa[0],"protocol_version") ) {
      send_gtp("= 2\n\n");
    } else if ( strstr(sa[0],"name")          ) {
      send_gtp("= argocorse_ichigo\n\n");
    } else if ( strstr(sa[0],"version")       ) {
      send_gtp("= 0.0.1\n\n");
    } else if ( strstr(sa[0],"list_commands" ) ) {
      send_gtp("= boardsize\nclear_board\nquit\nprotocol_version\nundo\n"
               "name\nversion\nlist_commands\nkomi\ngenmove\nplay\n\n");
    } else if ( strstr(sa[0],"komi") ) {
      komi = atof( sa[1] );
      send_gtp("= \n\n");
    } else if ( strstr(sa[0],"undo") ) {
      undo();
      send_gtp("= \n\n");
    } else if ( strstr(sa[0],"genmove") ) {
      int color = 1;
      if ( tolower(sa[1][0])=='w' ) color = 2;

      z = play_computer_move(color, SEARCH_UCT);
      send_gtp("= %s\n\n",get_char_z(z));
    } else if ( strstr(sa[0],"play") ) {  // "play b c4", "play w d17"
      int color = 1;
      if ( tolower(sa[1][0])=='w' ) color = 2;
      ax = tolower( sa[2][0] );
      x = ax - 'a' + 1;
      if ( ax >= 'i' ) x--;
      y = atoi(&sa[2][1]);
      z = get_z(x, B_SIZE-y+1);
      if ( strstr(sa[2],"pass") || strstr(sa[2],"PASS") ) z = 0;  // pass
      add_moves(z, color, 0);
      send_gtp("= \n\n");
    } else {
      send_gtp("? unknown_command\n\n");
    }
  }
}


int main()
{
	//srand( (unsigned)time( NULL ) );
	init_board();
	make_hashboard();

	if (0) { selfplay(); return 0; }
	if (0) { test_playout(); return 0; }

	gtp_loop();

//	getchar();
	return 0;
}
