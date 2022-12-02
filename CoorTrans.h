#include <cstdio>
#include <Eigen>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;



//�������ĵ��ĵع�����
struct coorECEF_withname
{
	double x;
	double y;
	double z;
	string name;
};

//���ĵع�����
struct coorECEF
{
	double x;
	double y;
	double z;
};

//�������
struct coorBLH
{
	double B;
	double L;
	double H;
};

//�������Ĵ������
struct coorBLH_withname
{
	double B;
	double L;
	double H;
	string name;
};

//վ������
struct coorNEU
{
	double N;
	double E;
	double U;
};

//��������վ������
struct coorNEU_withname
{
	double N;
	double E;
	double U;
	string name;
};

//�ڷ��Ͼ���ָ��
struct precisionPara
{
	int n;//����������
	int r;//�����
	double sigma0;//��λȨ�����
	vector<coorECEF_withname> e;
};

//����Ͼ���ָ��
struct precisionResidual
{
	double variance;//����
	vector<coorECEF_withname> e;
};

//�Ĳ���
struct parameter_4
{
	double dx;
	double dy;
	double R;
	double K;
	precisionPara precision;
};

//������
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

//�߲���
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

//ʮ������
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

//�����Ĳ���
bool computePara_4(map<string, vector<coorECEF>> coor, parameter_4 &para);

//��Դ����ת����Ŀ������ �Ĳ���
bool computeDstCoorPara_4(vector<coorECEF_withname>src, vector<coorECEF_withname>&dst,parameter_4 para);

//����������
bool computePara_6(map<string, vector<coorECEF>> coor, parameter_6& para);

//��Դ����ת����Ŀ������ ������
bool computeDstCoorPara_6(vector<coorECEF_withname>src, vector<coorECEF_withname>& dst, parameter_6 para);

//�����߲���
bool computePara_7(map<string, vector<coorECEF>> coor, parameter_7& para);

//��Դ����ת����Ŀ������ �߲���
bool computeDstCoorPara_7(vector<coorECEF_withname>src, vector<coorECEF_withname>& dst, parameter_7 para);

//����ʮ������
bool computePara_13(map<string, vector<coorECEF>> coor, parameter_13& para);

//��Դ����ת����Ŀ������ ʮ������
bool computeDstCoorPara_13(vector<coorECEF_withname>src, vector<coorECEF_withname>& dst, parameter_13 para);

//����в�
bool computeResidual(vector<coorECEF_withname>estimate, vector<coorECEF_withname> real, precisionResidual& res);

//XYZתBLH
bool coorXYZ2BLH(vector<coorECEF_withname>src, vector<coorBLH_withname>& dst);

//BLHתXYZ
bool coorBLH2XYZ(vector<coorBLH_withname> src, vector<coorECEF_withname>& dst);

//XYZתNEU
bool coorXYZ2NEU(coorECEF_withname center, vector<coorECEF_withname> src, vector<coorNEU_withname>& dst);