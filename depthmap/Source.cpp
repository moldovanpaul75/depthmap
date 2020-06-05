#include "stdafx.h"
#include "common.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <queue>
#include <random>

/*
*	Moldovan Paul Andrei
*	30235
*   Depth Map Estimation with Census Transform and Hamming Distance
*/

using namespace std;


uchar* censusTransform(Mat img, int h, int v) {
	
	imshow("original img", img);
	Mat imgTemp = Mat::zeros(img.rows, img.cols, CV_8UC1);	//imagine plina de 0

	unsigned int census = 0;
	int v1 = v / 2;
	int h1 = h / 2;	

	//calculez size la data
	uchar* data_;
	int pixelBitSize = v * h - 1;	//calcuzel cati biti o sa fie pe pixel, pt 3x3 => 8biti/pixel pentru ca ignoram pixelul din centru
	int totalPixels = img.rows*img.cols*pixelBitSize;	//numarul total de biti pentru toata imaginea
	int dataSize = (totalPixels + sizeof(uchar) - 1) / sizeof(uchar); //size-ul la data
	data_ = new uchar[dataSize]; //vector de uchar unde salvezi bitii
	memset(data_, 0, dataSize * sizeof(uchar)); //setez tot vectorul 0

	//printf("%d %d %d\n\n", pixelBitSize, totalPixels, dataSize);

	int count = 0;	
	int pixel_pos, bit_pos, bit_pos_start;
	uchar val, center_val;

	//parcurg imaginea pixel cu pixel incepand de la v1/h1 astfel ca sa nu am pading(sa nu iasa kernelul din imagine)
	for (int y = v1; y < img.rows - v1; y++) {
		for (int x = h1; x < img.cols - h1; x++) {

			census = 0;
			count = 0;

			pixel_pos = y * img.cols + x;
			bit_pos_start = pixel_pos * pixelBitSize; //pixel bit size

			center_val = img.at<uchar>(y, x);
			//parcurg kernelul
			for (int i = -v1; i <= v1; i++) {
				for (int j = -h1; j <= h1; j++) {

					if (i == 0 && j == 0) continue; //skip center

					int w_x = max(min(x + j, img.cols-1), 0);	//window x
					int w_y = max(min(y + i, img.rows-1), 0);	//window y

					census <<= 1;	//shift cu 1 la stanga (census = census << 1)
					bit_pos = bit_pos_start + count++;	//pentru pixelul curent
					val = img.at<uchar>(w_y, w_x);

					if (val < center_val) {	//verific daca intensitatea pixelului curent 
											//e mai mica decat cea a pixelului din centrul kernelului
						census += 0;
					}
					else {
						int shift = bit_pos % sizeof(uchar);				
						data_[bit_pos / sizeof(uchar)] |= (1 << shift);		
						//sau logic cu 1 shiftat pentru pozitia curenta

						census += 1;
					}
				}
			}
			//printf("%d %d %d\n", bit_pos, data_[bit_pos / sizeof(uchar)], census);

			imgTemp.at<uchar>(y, x) = census;
		}
	}
	imshow("census", imgTemp);
	return data_;
}

int getBitVal(uchar *data, int pos) {
	uchar elem = data[pos / sizeof(uchar)];
	int shift = pos % sizeof(uchar);
	elem = (elem >> shift);
	return elem & 1;
}


int getDistance(uchar *data_, int y1, int x1, int y2, int x2, int h, int v, int cols, int rows){
	int v1 = v / 2;
	int h1 = h / 2;
	int pixelBitSize = v * h - 1;

	int bit_pos_start1 = (y1 * cols + x1) * pixelBitSize;
	int bit_pos_start2 = (y2 * cols + x2) * pixelBitSize;

	int dis = 0;
	int count = 0;
	for (int i = -v1; i <= v1; i++){
		for (int j = -h1; j <= h1; j++){

			int bit_pos1 = bit_pos_start1 + count;
			int bit_pos2 = bit_pos_start2 + count;

			count++;
			int bit_val1 = getBitVal(data_, bit_pos1);
			int bit_val2 = getBitVal(data_, bit_pos2);
			printf("%d %d %d %d\n", bit_pos1, bit_pos2, bit_val1, bit_val2);
			dis += bit_val1 != bit_val2;	//pt fiecare bit diferit distanta++
		}
	}
	return dis;
}

int getDistance2(uchar *data_left, uchar *data_right, int y1, int x1, int y2, int x2, int h, int v, int cols, int rows) {
	int v1 = v / 2;
	int h1 = h / 2;
	int pixelBitSize = v * h - 1;

	int bit_pos_start1 = (y1 * cols + x1) * pixelBitSize;
	int bit_pos_start2 = (y2 * cols + x2) * pixelBitSize;

	int dis = 0;
	int count = 0;
	for (int i = -v1; i <= v1; i++) {
		for (int j = -h1; j <= h1; j++) {

			int bit_pos1 = bit_pos_start1 + count;
			int bit_pos2 = bit_pos_start2 + count;

			count++;
			int bit_val1 = getBitVal(data_left, bit_pos1);
			int bit_val2 = getBitVal(data_right, bit_pos2);
			//printf("%d %d %d %d\n", bit_pos1, bit_pos2, bit_val1, bit_val2);
			dis += bit_val1 != bit_val2;
		}
	}
	return dis;
}



void computeDispMap(Mat left_img, Mat right_img, int min_disparity, int max_disparity, int w) {

	int k = w / 2;
	Mat disp_map = Mat(left_img.rows, left_img.cols, CV_8UC1);	//rezultatul

	uchar *data_left = censusTransform(left_img, w, w);
	uchar *data_right = censusTransform(right_img, w, w);

	int count = 0, dis = 0;
	float cost = 0;
	long min_cost = 0;
	
	//parcug imaginile
	for (int y = 0; y < left_img.rows; y++) {
		for (int x = 0; x < left_img.cols; x++) {

			min_cost = 1e9; //10^9
			dis = 0;
			//pentru fiecare disparitate(diferenta)
			for (int d = min_disparity; d <= max_disparity; d++) {
				
				//daca nu e in imagine continue
				if (x - d < 0) continue;

				cost = 0;
				count = 0;
				for (int i = -k; i <= k; i++) {			//parcurg kernelul
					for (int j = -k; j <= k; j++) {
						
						int w_x = x + j;				//window x
						int w_y = y + i;				//window y
						int w_x2 = w_x - d;	//current x - disparity -> img2 x

						//daca nu e in imagine
						if (w_x < 0 || w_y < 0 || w_x >= left_img.cols || w_y >= left_img.rows || w_x2 < 0 || w_x2 >= left_img.cols) continue;

						//calculez costul cu distanta hamming
						cost += getDistance2(data_left, data_right, w_y, w_x, w_y, w_x2, w, w, left_img.cols, left_img.rows);
						count += 1;

					}
				}

				cost = cost/count;
				if (cost < min_cost) {
					min_cost = cost;
					dis = d * 3;
				}
			}
			disp_map.at<uchar>(y, x) = dis;
		}
	}

	imshow("disparity", disp_map);
	waitKey(0);
}

void computeDispMap2(Mat left_img, Mat right_img, int ndisp, int w, int w2) {

	int k = w / 2;
	Mat disp_map = Mat(left_img.rows, left_img.cols, CV_8UC1);

	//calculez transformatele pentru cele doua imagini
	uchar *data_left = censusTransform(left_img, w2, w2);		
	uchar *data_right = censusTransform(right_img, w2, w2);

	int dis = 0;
	float cost = 0;
	long min_cost = 0;
	int w_x, w_y, w_x2;
	//parcug imaginile
	for (int y = ndisp; y < left_img.rows  - ndisp; y++) {
		for (int x = ndisp; x < left_img.cols - ndisp; x++) {

			min_cost = 1e9; //10^9, asignez o valoare foarte mare pentru costul minim
			dis = 0;
			//pentru fiecare disparitate(diferenta)
			for (int d = x-ndisp; d <= x; d++) {	//ndisp = disparitatea maxima

				cost = 0;
				for (int i = -w; i <= w; i++) {			//parcurg kernelul
					for (int j = -w; j <= w; j++) {

						w_x = x + j;				//x img1
						w_y = y + i;				//y img1
						w_x2 = j + d;				//x2 img2

						//daca nu e in imagine
						if (w_x < 0 || w_y < 0 || w_x >= left_img.cols || w_y >= left_img.rows || w_x2 < 0 || w_x2 >= left_img.cols) continue;

						//calculez distanta Hamming dintre cele doua puncte
						cost += getDistance2(data_left, data_right, w_y, w_x, w_y, w_x2, w2, w2, left_img.cols, left_img.rows);

						//folosind distanta Manhattan
						//cost += abs(left_img.at<uchar>(w_y, w_x) - right_img.at<uchar>(w_y, w_x2));
					}
				}
				//daca costul calculat e < min_cost atunci devine min cost si avem o noua disparitate d
				if (cost < min_cost) {
					min_cost = cost;
					dis = d;
				}
			}
			disp_map.at<uchar>(y, x) = 3*abs(x - dis) * (255. / ndisp);	//formula pentru inversul disparitatii (*3 pentru ca imaginea sa fie mai deschisa)
		}
	}

	imshow("disparity", disp_map);
	waitKey(0);
}



int getPoint(int x, int y, Mat left_img, Mat right_img, int ndisp, int w, uchar *data_left, uchar *data_right) {
	long min_cost = 1e9;  //10^9
	int chosen_x = 0;

	for (int k = x - ndisp; k <= x; k++) {
		long error = 0;
		for (int i = -w; i <= w; i++) {
			for (int j = -w; j <= w; j++) {

				int w_x = x + j;
				int w_y = y + i;
				int w_x2 = k + j;
				
				if (w_x < 0 || w_y < 0 || w_x >= left_img.cols || w_y >= left_img.rows || w_x2 < 0 || w_x2 >= left_img.cols) continue;
				error += abs(left_img.at<uchar>(w_y, w_x) - right_img.at<uchar>(w_y, w_x2));
			}
		}
		if (error < min_cost) {
			min_cost = error;
			chosen_x = k;
		}

	}
	return chosen_x;
}

void computeDisparityMap(Mat left_img, Mat right_img, int ndisp, int w) {

	int k = w / 2;
	Mat disp_map = Mat(left_img.rows, left_img.cols, CV_8UC1);	//rezultatul

	uchar *data_left = censusTransform(left_img, 3, 3);
	uchar *data_right = censusTransform(right_img, 3, 3);

	for (int y = ndisp; y < left_img.rows - ndisp; y++) {
		for (int x = ndisp; x < left_img.cols - ndisp; x++) {

			int right_x = getPoint(x, y, left_img, right_img, ndisp, w, data_left, data_right);
			int dis = abs(x - right_x);
			disp_map.at<uchar>(y, x) = 3 * dis *  (255. / ndisp);
		}
	}
	imshow("disparity", disp_map);
}



void hamming_distance_demo() {
	int datas[6][5] = {
		{ 0, 0, 0, 0, 0 },
		{ 0, 0, 1, 0, 0 },
		{ 0, 0, 0, 0, 0 },
		{ 0, 1, 1, 1, 0 },
		{ 0, 1, 0, 1, 0 },
		{ 0, 1, 1, 1, 0 },
	};
	printf("demo distanta hamming\n");
	Mat image1 = cv::Mat::ones(6, 5, CV_8UC1);
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			image1.at<uchar>(i, j) = datas[i][j];
		}
	}
	uchar *data = censusTransform(image1, 3, 3);
	int dis = getDistance(data, 1, 2, 4, 2, 3, 3, 5, 5);
	printf("rezultat: %d \n", dis);
}


int main() {

	Mat left_img = imread("Images/pentagon_left.bmp", 0);
	Mat right_img = imread("Images/pentagon_right.bmp", 0);


	//hamming_distance_demo();
	//censusTransform(left_img, 7, 7);
	//computeDispMap2(left_img, right_img, 20, 3, 3);			//hamming distance test
	computeDisparityMap(left_img, right_img, 20, 7);		//manhattan distance test


	//getchar();
	waitKey(0);
	return 0;
}