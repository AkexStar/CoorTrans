#include <cstdio>
#include <Eigen>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;



//带点名的地心地固坐标
struct coorCECF_withname
{
	double x;
	double y;
	double z;
	string name;
};

//地心地固坐标
struct coorCECF
{
	double x;
	double y;
	double z;
};

//大地坐标
struct coorBLH
{
	double B;
	double L;
	double H;
};

//带点名的大地坐标
struct coorBLH_withname
{
	double B;
	double L;
	double H;
	string name;
};

//站心坐标
struct coorNEU
{
	double N;
	double E;
	double U;
};

//带点名的站心坐标
struct coorNEU_withname
{
	double N;
	double E;
	double U;
	string name;
};

//内符合精度指标
struct precisionPara
{
	int n;//条件方程数
	int r;//冗余度
	double sigma0;//单位权中误差
	vector<coorCECF_withname> e;
};

//外符合精度指标
struct precisionResidual
{
	double variance;//方差
	vector<coorCECF_withname> e;
};

//四参数
struct parameter_4
{
	double dx;
	double dy;
	double R;
	double K;
	precisionPara precision;
};

//六参数
struct parameter_6
{
	double dx;
	double dy;
	double dz;
	double Rx;
	double Ry;
	double Rz;
	precisionPara precision;
};

//七参数
struct parameter_7
{
	double dx;
	double dy;
	double dz;
	double Rx;
	double Ry;
	double Rz;
	double K;
	precisionPara precision;
};

//十三参数
struct parameter_13
{
	double dx;
	double dy;
	double dz;
	Matrix<double, 3, 3> R;
	double K;
	int count;
	precisionPara precision;
};

//计算四参数
bool computePara_4(map<string, vector<coorCECF>> coor, parameter_4 &para);

//将源坐标转换到目标坐标 四参数
bool computeDstCoorPara_4(vector<coorCECF_withname>src, vector<coorCECF_withname>&dst,parameter_4 para);

//计算六参数
bool computePara_6(map<string, vector<coorCECF>> coor, parameter_6& para);

//将源坐标转换到目标坐标 六参数
bool computeDstCoorPara_6(vector<coorCECF_withname>src, vector<coorCECF_withname>& dst, parameter_6 para);

//计算七参数
bool computePara_7(map<string, vector<coorCECF>> coor, parameter_7& para);

//将源坐标转换到目标坐标 七参数
bool computeDstCoorPara_7(vector<coorCECF_withname>src, vector<coorCECF_withname>& dst, parameter_7 para);

//计算十三参数
bool computePara_13(map<string, vector<coorCECF>> coor, parameter_13& para);

//将源坐标转换到目标坐标 十三参数
bool computeDstCoorPara_13(vector<coorCECF_withname>src, vector<coorCECF_withname>& dst, parameter_13 para);

//计算残差
bool computeResidual(vector<coorCECF_withname>estimate, vector<coorCECF_withname> real, precisionResidual& res);

//XYZ转BLH
bool coorXYZ2BLH(vector<coorCECF_withname>src, vector<coorBLH_withname>& dst);

//BLH转XYZ
bool coorBLH2XYZ(vector<coorBLH_withname> src, vector<coorCECF_withname>& dst);

//XYZ转NEU
bool coorXYZ2NEU(coorCECF_withname center, vector<coorCECF_withname> src, vector<coorNEU_withname>& dst);