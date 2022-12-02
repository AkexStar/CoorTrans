#include <cstdio>
#include <Eigen>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;



//�������ĵ��ĵع�����
struct coorCECF_withname
{
	double x;
	double y;
	double z;
	string name;
};

//���ĵع�����
struct coorCECF
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
	vector<coorCECF_withname> e;
};

//����Ͼ���ָ��
struct precisionResidual
{
	double variance;//����
	vector<coorCECF_withname> e;
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
bool computePara_4(map<string, vector<coorCECF>> coor, parameter_4 &para);

//��Դ����ת����Ŀ������ �Ĳ���
bool computeDstCoorPara_4(vector<coorCECF_withname>src, vector<coorCECF_withname>&dst,parameter_4 para);

//����������
bool computePara_6(map<string, vector<coorCECF>> coor, parameter_6& para);

//��Դ����ת����Ŀ������ ������
bool computeDstCoorPara_6(vector<coorCECF_withname>src, vector<coorCECF_withname>& dst, parameter_6 para);

//�����߲���
bool computePara_7(map<string, vector<coorCECF>> coor, parameter_7& para);

//��Դ����ת����Ŀ������ �߲���
bool computeDstCoorPara_7(vector<coorCECF_withname>src, vector<coorCECF_withname>& dst, parameter_7 para);

//����ʮ������
bool computePara_13(map<string, vector<coorCECF>> coor, parameter_13& para);

//��Դ����ת����Ŀ������ ʮ������
bool computeDstCoorPara_13(vector<coorCECF_withname>src, vector<coorCECF_withname>& dst, parameter_13 para);

//����в�
bool computeResidual(vector<coorCECF_withname>estimate, vector<coorCECF_withname> real, precisionResidual& res);

//XYZתBLH
bool coorXYZ2BLH(vector<coorCECF_withname>src, vector<coorBLH_withname>& dst);

//BLHתXYZ
bool coorBLH2XYZ(vector<coorBLH_withname> src, vector<coorCECF_withname>& dst);

//XYZתNEU
bool coorXYZ2NEU(coorCECF_withname center, vector<coorCECF_withname> src, vector<coorNEU_withname>& dst);