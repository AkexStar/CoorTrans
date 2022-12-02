#include "CoorTrans.h"

//计算四参数
bool computePara_4(map<string, vector<coorECEF>> coor, parameter_4 &para)
{
	vector<coorECEF> src;	//存储源坐标系坐标
	vector<coorECEF> dst;	//存储目标坐标系坐标
	vector<string> name;	//储存点名
	for (auto it = coor.begin(); it != coor.end(); ++it)
	{
		if(it->second.size()==2)
		{
			src.push_back(it->second[0]);
			dst.push_back(it->second[1]);
			name.push_back(it->first);
		}
	}
		
	int src_n = src.size();
	int dst_n = dst.size();
	if(src_n!=dst_n)
	{
		cout << "源坐标系点数与目标坐标系点数不相等！\n" << endl;
		return false;
	}
	if(src_n*2 - 4<0)
	{
		cout << "条件方程个数不足！无法求解！\n";
		return false;
	}
	int d = src_n;//点数量
	int n=para.precision.n = src_n * 2;//条件方程数量
	para.precision.r = para.precision.n - 4;//冗余度
	
	MatrixXd B(n,4);//系数矩阵
	Matrix<double, Dynamic, 1> l;//常数向量
	l.resize(n, 1);
	MatrixXd x_s(d, 1), y_s(d, 1), x_t(d, 1), y_t(d, 1);
	for (int i = 0; i < d; i++)
	{
		x_s(i, 0) = src[i].x;
		y_s(i, 0) = src[i].y;
		x_t(i, 0) = dst[i].x;
		y_t(i, 0) = dst[i].y;
	}
	l << x_t - x_s,
		 y_t - y_s;
	
	B << MatrixXd::Ones(d, 1), MatrixXd::Zero(d, 1), y_s, x_s,
		MatrixXd::Zero(d, 1), MatrixXd::Ones(d, 1), -x_s, y_s;
	Matrix<double, 4, 1> x;
	x = (B.transpose() * B).inverse() * B.transpose() * l;//求解法方程
	para.dx = x(0, 0);
	para.dy = x(1, 0);
	para.R = x(2, 0);
	para.K = x(3, 0);
	MatrixXd V(n,1);
	//计算改正数V
	V = - (MatrixXd::Identity(n,n)- B * (B.transpose() * B).inverse() * B.transpose())*l;
	for (int i=0;i<d;i++)
	{
		coorECEF_withname temp_e={V(i,0),V(i+d,0),0,name[i]};
		para.precision.e.push_back(temp_e);
	}
	//计算sigma0
	para.precision.sigma0 = ((V.transpose() * V) / para.precision.r)(0,0);
	para.precision.sigma0 = sqrt(para.precision.sigma0);
	return true;
}

bool computeDstCoorPara_4(vector<coorECEF_withname>src, vector<coorECEF_withname>& dst, parameter_4 para)
{
	MatrixXd x(4, 1);
	x(0, 0) = para.dx;
	x(1, 0) = para.dy;
	x(2, 0) = para.R;
	x(3, 0) = para.K;
	int d = src.size();
	if(d==0)
	{
		cout << "源坐标数据文件为空！\n";
		return false;
	}
	int n = 2 * d;
	MatrixXd B(n, 4);
	MatrixXd x_s(d, 1), y_s(d, 1), x_t(d, 1), y_t(d, 1);
	for (int i = 0; i < d; i++)
	{
		x_s(i, 0) = src[i].x;
		y_s(i, 0) = src[i].y;
	}
	B << MatrixXd::Ones(d, 1), MatrixXd::Zero(d, 1), y_s, x_s,
		MatrixXd::Zero(d, 1), MatrixXd::Ones(d, 1), -x_s, y_s;
	MatrixXd y(n, 1), temp(n,1);
	temp << x_s,
		y_s;
	y = temp + B * x;
	for(int i=0;i<d;i++)
	{
		coorECEF_withname temp={y(i,0),y(i+d,0),0,src[i].name};
		dst.push_back(temp);
	}
	return true;
}

//计算七参数
bool computePara_7(map<string, vector<coorECEF>> coor, parameter_7& para)
{
	vector<coorECEF> src;	//存储源坐标系坐标
	vector<coorECEF> dst;	//存储目标坐标系坐标
	vector<string> name;	//储存点名
	for (auto it = coor.begin(); it != coor.end(); ++it)
	{
		if (it->second.size() == 2)
		{
			src.push_back(it->second[0]);
			dst.push_back(it->second[1]);
			name.push_back(it->first);
		}
	}

	int src_n = src.size();
	int dst_n = dst.size();
	if (src_n != dst_n)
	{
		cout << "源坐标系点数与目标坐标系点数不相等！\n" << endl;
		return false;
	}
	if (src_n * 3 - 6 < 0)
	{
		cout << "条件方程个数不足！无法求解！\n";
		return false;
	}
	int d = src_n;//点数量
	int n = para.precision.n = src_n * 3;//条件方程数量
	para.precision.r = para.precision.n - 7;//冗余度

	MatrixXd B(n, 7);//系数矩阵
	MatrixXd l(n, 1);//常数向量
	MatrixXd x_s(d, 1), y_s(d, 1), z_s(d,1),
	x_t(d, 1), y_t(d, 1), z_t(d,1);
	for (int i = 0; i < d; i++)
	{
		x_s(i, 0) = src[i].x;
		y_s(i, 0) = src[i].y;
		z_s(i, 0) = src[i].z;
		x_t(i, 0) = dst[i].x;
		y_t(i, 0) = dst[i].y;
		z_t(i, 0) = dst[i].z;
	}
	l << x_t - x_s,
		y_t - y_s,
		z_t - z_s;
	MatrixXd ONE = MatrixXd::Ones(d, 1);
	MatrixXd ZERO = MatrixXd::Zero(d, 1);
	B << ONE, ZERO, ZERO, ZERO, -z_s, y_s, x_s,
		ZERO, ONE, ZERO, z_s, ZERO, -x_s, y_s,
		ZERO, ZERO, ONE, -y_s, x_s, ZERO, z_s;
	Matrix<double, 7, 1> x;
	x = (B.transpose() * B).inverse() * B.transpose() * l;//求解法方程
	para.dx = x(0, 0);
	para.dy = x(1, 0);
	para.dz = x(2, 0);
	para.Rx = x(3, 0);
	para.Ry = x(4, 0);
	para.Rz = x(5, 0);
	para.K = x(6, 0);
	MatrixXd V(n, 1);
	//计算改正数V
	V = - (MatrixXd::Identity(n, n) - B * (B.transpose() * B).inverse() * B.transpose()) * l;
	for (int i = 0; i < d; i++)
	{
		coorECEF_withname temp_e = { V(i,0),V(i + d,0), V(i+2*d,0),name[i] };
		para.precision.e.push_back(temp_e);
	}
	//计算sigma0
	para.precision.sigma0 = ((V.transpose() * V) / para.precision.r)(0, 0);
	para.precision.sigma0 = sqrt(para.precision.sigma0);
	return true;
}

bool computeDstCoorPara_7(vector<coorECEF_withname>src, vector<coorECEF_withname>& dst, parameter_7 para)
{
	MatrixXd x(7, 1);
	x(0, 0) = para.dx;
	x(1, 0) = para.dy;
	x(2, 0) = para.dz;
	x(3, 0) = para.Rx;
	x(4, 0) = para.Ry;
	x(5, 0) = para.Rz;
	x(6, 0) = para.K;
	int d = src.size();
	if (d == 0)
	{
		cout << "源坐标数据文件为空！\n";
		return false;
	}
	int n = 3 * d;
	MatrixXd B(n, 7);
	MatrixXd x_s(d, 1), y_s(d, 1), z_s(d, 1);
	MatrixXd x_t(d, 1), y_t(d, 1), z_t(d,1);
	for (int i = 0; i < d; i++)
	{
		x_s(i, 0) = src[i].x;
		y_s(i, 0) = src[i].y;
		z_s(i, 0) = src[i].z;
	}
	MatrixXd ONE = MatrixXd::Ones(d, 1);
	MatrixXd ZERO = MatrixXd::Zero(d, 1);
	B << ONE, ZERO, ZERO, ZERO, -z_s, y_s, x_s,
		ZERO, ONE, ZERO, z_s, ZERO, -x_s, y_s,
		ZERO, ZERO, ONE, -y_s, x_s, ZERO, z_s;
	
	MatrixXd y(n, 1), temp(n, 1);
	temp << x_s,
		y_s,
		z_s;
	y = temp + B * x;
	for (int i = 0; i < d; i++)
	{
		coorECEF_withname temp = { y(i,0),y(i + d,0),y(i+2*d,0),src[i].name };
		dst.push_back(temp);
	}
	return true;
}

//计算六参数
bool computePara_6(map<string, vector<coorECEF>> coor, parameter_6& para)
{
	vector<coorECEF> src;	//存储源坐标系坐标
	vector<coorECEF> dst;	//存储目标坐标系坐标
	vector<string> name;	//储存点名
	for (auto it = coor.begin(); it != coor.end(); ++it)
	{
		if (it->second.size() == 2)
		{
			src.push_back(it->second[0]);
			dst.push_back(it->second[1]);
			name.push_back(it->first);
		}
	}

	int src_n = src.size();
	int dst_n = dst.size();
	if (src_n != dst_n)
	{
		cout << "源坐标系点数与目标坐标系点数不相等！\n" << endl;
		return false;
	}
	if (src_n * 3 - 6 < 0)
	{
		cout << "条件方程个数不足！无法求解！\n";
		return false;
	}
	int d = src_n;//点数量
	int n = para.precision.n = src_n * 3;//条件方程数量
	para.precision.r = para.precision.n - 6;//冗余度

	MatrixXd B(n, 6);//系数矩阵
	MatrixXd l(n, 1);//常数向量
	MatrixXd x_s(d, 1), y_s(d, 1), z_s(d, 1),
		x_t(d, 1), y_t(d, 1), z_t(d, 1);
	for (int i = 0; i < d; i++)
	{
		x_s(i, 0) = src[i].x;
		y_s(i, 0) = src[i].y;
		z_s(i, 0) = src[i].z;
		x_t(i, 0) = dst[i].x;
		y_t(i, 0) = dst[i].y;
		z_t(i, 0) = dst[i].z;
	}
	l << x_t - x_s,
		y_t - y_s,
		z_t - z_s;
	MatrixXd ONE = MatrixXd::Ones(d, 1);
	MatrixXd ZERO = MatrixXd::Zero(d, 1);
	B << ONE, ZERO, ZERO, ZERO, -z_s, y_s,
		ZERO, ONE, ZERO, z_s, ZERO, -x_s,
		ZERO, ZERO, ONE, -y_s, x_s, ZERO;
	Matrix<double, 6, 1> x;
	x = (B.transpose() * B).inverse() * B.transpose() * l;//求解法方程
	para.dx = x(0, 0);
	para.dy = x(1, 0);
	para.dz = x(2, 0);
	para.Rx = x(3, 0);
	para.Ry = x(4, 0);
	para.Rz = x(5, 0);
	MatrixXd V(n, 1);
	//计算改正数V
	V = -(MatrixXd::Identity(n, n) - B * (B.transpose() * B).inverse() * B.transpose()) * l;
	for (int i = 0; i < d; i++)
	{
		coorECEF_withname temp_e = { V(i,0),V(i + d,0), V(i + 2 * d,0),name[i] };
		para.precision.e.push_back(temp_e);
	}
	//计算sigma0
	para.precision.sigma0 = ((V.transpose() * V) / para.precision.r)(0, 0);
	para.precision.sigma0 = sqrt(para.precision.sigma0);
	return true;
}

bool computeDstCoorPara_6(vector<coorECEF_withname>src, vector<coorECEF_withname>& dst, parameter_6 para)
{
	MatrixXd x(6, 1);
	x(0, 0) = para.dx;
	x(1, 0) = para.dy;
	x(2, 0) = para.dz;
	x(3, 0) = para.Rx;
	x(4, 0) = para.Ry;
	x(5, 0) = para.Rz;
	int d = src.size();
	if (d == 0)
	{
		cout << "源坐标数据文件为空！\n";
		return false;
	}
	int n = 3 * d;
	MatrixXd B(n, 6);
	MatrixXd x_s(d, 1), y_s(d, 1), z_s(d, 1);
	MatrixXd x_t(d, 1), y_t(d, 1), z_t(d, 1);
	for (int i = 0; i < d; i++)
	{
		x_s(i, 0) = src[i].x;
		y_s(i, 0) = src[i].y;
		z_s(i, 0) = src[i].z;
	}
	MatrixXd ONE = MatrixXd::Ones(d, 1);
	MatrixXd ZERO = MatrixXd::Zero(d, 1);
	B << ONE, ZERO, ZERO, ZERO, -z_s, y_s,
		ZERO, ONE, ZERO, z_s, ZERO, -x_s,
		ZERO, ZERO, ONE, -y_s, x_s, ZERO;

	MatrixXd y(n, 1), temp(n, 1);
	temp << x_s,
		y_s,
		z_s;
	y = temp + B * x;
	for (int i = 0; i < d; i++)
	{
		coorECEF_withname temp = { y(i,0),y(i + d,0),y(i + 2 * d,0),src[i].name };
		dst.push_back(temp);
	}
	return true;
}

//克罗内克积
extern MatrixXd Kron(Eigen::MatrixXd m1, Eigen::MatrixXd m2) {
	int m1R, m1C, m2R, m2C;
	m1R = m1.rows();
	m1C = m1.cols();

	m2R = m2.rows();
	m2C = m2.cols();

	MatrixXd m3(m1R * m2R, m1C * m2C);

	for (int i = 0; i < m1C; i++) {
		for (int j = 0; j < m1R; j++) {
			m3.block(i * m2R, j * m2C, m2R, m2C) = m1(i, j) * m2;
		}
	}
	return m3;
}

//计算十三参数
//bool computePara_13(map<string, vector<coorECEF>> coor, parameter_13& para)
//{
//	vector<coorECEF> src;	//存储源坐标系坐标
//	vector<coorECEF> dst;	//存储目标坐标系坐标
//	vector<string> name;	//储存点名
//	for (auto it = coor.begin(); it != coor.end(); ++it)
//	{
//		if (it->second.size() == 2)
//		{
//			src.push_back(it->second[0]);
//			dst.push_back(it->second[1]);
//			name.push_back(it->first);
//		}
//	}
//
//	int src_n = src.size();
//	int dst_n = dst.size();
//	if (src_n != dst_n)
//	{
//		cout << "源坐标系点数与目标坐标系点数不相等！\n" << endl;
//		return false;
//	}
//	if (src_n * 3 - 7  < 0)
//	{
//		cout << "条件方程个数不足！无法求解！\n";
//		return false;
//	}
//	
//	int d = src_n;//点数量
//	int n = para.precision.n = src_n * 3;//观测方程数量
//	para.precision.r = para.precision.n - 7;//冗余度
//	int s = 5;//限制条件数
//	int u = 12;//参数个数
//	
//	MatrixXd x_s(d, 1), y_s(d, 1), z_s(d, 1),
//		x_t(d, 1), y_t(d, 1), z_t(d, 1);
//	for (int i = 0; i < d; i++)
//	{
//		x_s(i, 0) = src[i].x;
//		y_s(i, 0) = src[i].y;
//		z_s(i, 0) = src[i].z;
//		x_t(i, 0) = dst[i].x;
//		y_t(i, 0) = dst[i].y;
//		z_t(i, 0) = dst[i].z;
//	}
//
//	MatrixXd y(n, 1);
//	y << x_t,//观测向量
//		y_t,
//		z_t;
//	MatrixXd A(n, u);//设计矩阵
//	MatrixXd P = MatrixXd::Identity(n, n);//权阵
//	MatrixXd temp(d, 3);
//	temp << x_s, y_s, z_s;
//	MatrixXd EYE3 = MatrixXd::Identity(3,3);
//	MatrixXd ONEd = MatrixXd::Ones(d, 1);
//	A << Kron(EYE3, temp),Kron(EYE3, ONEd);
//	
//	MatrixXd ksai0(u, 1),ksai1(u,1);
//	ksai0 = (A.transpose() * A).inverse() * A.transpose() * y;
//	ksai1 = ksai0 + MatrixXd::Ones(u, 1);
//	MatrixXd e(n, 1);
//	int count = 0;
//	while((ksai1-ksai0).lpNorm<2>()>1e-10)
//	{
//		count++;
//		ksai1 = ksai0;
//		MatrixXd ksai = ksai0;
//		double ksai11 = ksai(0, 0);
//		double ksai12 = ksai(1, 0);
//		double ksai13 = ksai(2, 0);
//		double ksai21 = ksai(3, 0);
//		double ksai22 = ksai(4, 0);
//		double ksai23 = ksai(5, 0);
//		double ksai31 = ksai(6, 0);
//		double ksai32 = ksai(7, 0);
//		double ksai33 = ksai(8, 0);
//
//		MatrixXd C(s,u);
//		C<< 2 * ksai11, 0, -2 * ksai31, 2 * ksai12, 0, -2 * ksai32, 2 * ksai13, 0, -2 * ksai33, 0, 0, 0,
//			0, 2 * ksai21, -2 * ksai31, 0, 2 * ksai22, -2 * ksai32, 0, 2 * ksai23, -2 * ksai33, 0, 0, 0,
//			ksai21, ksai11, 0, ksai22, ksai12, 0, ksai23, ksai13, 0, 0, 0, 0,
//			ksai31, 0, ksai11, ksai32, 0, ksai12, ksai33, 0, ksai13, 0, 0, 0,
//			0, ksai13, ksai21, 0, ksai32, ksai22, 0, ksai33, ksai23, 0, 0, 0;
//
//		MatrixXd c(s, 1);
//		c << ksai11 * ksai11 + ksai12 * ksai12 + ksai13 * ksai13 - ksai31 * ksai31 - ksai32 * ksai32 - ksai33 * ksai33,
//			ksai21* ksai21 + ksai22 * ksai22 + ksai23 * ksai23 - ksai31 * ksai31 - ksai32 * ksai32 - ksai33 * ksai33,
//			ksai11* ksai21 + ksai12 * ksai22 + ksai13 * ksai23,
//			ksai11* ksai31 + ksai12 * ksai32 + ksai13 * ksai33,
//			ksai31* ksai21 + ksai32 * ksai22 + ksai33 * ksai23;
//
//		MatrixXd LEFT(u + s, u + s);
//		LEFT<< (A.transpose()*P * A ), C.transpose(),
//			C, MatrixXd::Zero(s,s);
//			
//		MatrixXd RIGHT(u + s, 1);
//		RIGHT<< (A.transpose()*P* y),
//				C* ksai0 - c;
//		
//		MatrixXd x = LEFT.inverse() * RIGHT;
//		ksai0 = x.topRows(u);
//	}
//	para.count = count;
//	double ksai11 = ksai0(0, 0);
//	double ksai12 = ksai0(1, 0);
//	double ksai13 = ksai0(2, 0);
//	double ksai21 = ksai0(3, 0);
//	double ksai22 = ksai0(4, 0);
//	double ksai23 = ksai0(5, 0);
//	double ksai31 = ksai0(6, 0);
//	double ksai32 = ksai0(7, 0);
//	double ksai33 = ksai0(8, 0);
//	
//	para.K = sqrt(1/3.0*(ksai11* ksai11 + ksai12 * ksai12 + ksai13 * ksai13+
//		ksai21 * ksai21 + ksai22 * ksai22 + ksai23 * ksai23+
//		ksai31 * ksai31 + ksai32 * ksai32 + ksai33 * ksai33));
//	cout << e<< endl;
//	e = A * ksai0 - y;
//	cout << e<< endl;
//	para.precision.sigma0 = (e.transpose() * P * e)(0, 0) / para.precision.r;
//	para.precision.sigma0 = sqrt(para.precision.sigma0);
//	cout << ksai0 << endl;
//	para.R << ksai0.topRows(3).transpose(),
//			ksai0.middleRows(4, 3).transpose(),
//			ksai0.middleRows(7, 3).transpose();
//	para.R = para.R / para.K;
//	para.dx = ksai0(9, 0);
//	para.dy = ksai0(10, 0);
//	para.dz = ksai0(11, 0);
//	for (int i = 0; i < d; i++)
//	{
//		coorECEF_withname temp_e = { e(i,0),e(i + d,0), e(i + 2 * d,0),name[i] };
//		para.precision.e.push_back(temp_e);
//	}
//	return true;
//}
//

//计算十三参数
bool computePara_13(map<string, vector<coorECEF>> coor, parameter_13& para)
{
	vector<coorECEF> src;	//存储源坐标系坐标
	vector<coorECEF> dst;	//存储目标坐标系坐标
	vector<string> name;	//储存点名
	for (auto it = coor.begin(); it != coor.end(); ++it)
	{
		if (it->second.size() == 2)
		{
			src.push_back(it->second[0]);
			dst.push_back(it->second[1]);
			name.push_back(it->first);
		}
	}
	int src_d = src.size();//源坐标点数
	int dst_d = dst.size();//目标坐标点数
	if (src_d != dst_d)
	{
		cout << "源坐标系点数与目标坐标系点数不相等！\n" << endl;
		return false;
	}
	int d = src_d;//点数量
	int n = para.precision.n = src_d * 3;//条件方程数量3d
	int u = 13;//参数个数13
	int s = 6;//限制方程个数6
	
	if (src_d * 3 + s - u < 0)
	{
		cout << "条件方程个数不足！无法求解！\n";
		return false;
	}
	
	para.precision.r = para.precision.n + s -u;//冗余度
	
	MatrixXd x_s(d, 1), y_s(d, 1), z_s(d, 1),
		x_t(d, 1), y_t(d, 1), z_t(d, 1);//坐标向量
	for (int i = 0; i < d; i++)
	{
		x_s(i, 0) = src[i].x;
		y_s(i, 0) = src[i].y;
		z_s(i, 0) = src[i].z;
		x_t(i, 0) = dst[i].x;
		y_t(i, 0) = dst[i].y;
		z_t(i, 0) = dst[i].z;
	}

	MatrixXd P(n + s, n + s);//权阵 n+s * n+s
	P = MatrixXd::Identity(n + s, n + s);
	MatrixXd A(n, u);//误差方程系数 n*u
	MatrixXd l(n, 1);//误差方程常数向量 n*1
	MatrixXd x0(u, 1);//参数改正数初值向量 u*1
	MatrixXd x1(u, 1);//参数改正数向量 u*1
	MatrixXd EYE3 = MatrixXd::Identity(3, 3);//单位阵I3
	MatrixXd ONEd = MatrixXd::Ones(d, 1);//一列1 d行
	MatrixXd ZERO = MatrixXd::Zero(d, 1);//一列0 d行
	MatrixXd xyz_s(d, 3);//源坐标 并排x y z
	MatrixXd xyz_t(n, 1);//目标坐标 竖排x\n y\n z\n
	double a1, a2, a3, b1, b2, b3, c1, c2, c3, K, dx, dy, dz;//参数值
	MatrixXd R0(3, 3);//旋转矩阵初值
	xyz_s << x_s, y_s, z_s;
	xyz_t << x_t,
		y_t,
		z_t;
	MatrixXd X(u, 1);//参数值
	//*****计算参数的初值*****//
	{
		A.setZero(n, u - 1);
		A << Kron(EYE3, xyz_s), Kron(EYE3, ONEd);
		MatrixXd X0(12, 1);
		X0 = (A.transpose() * A).inverse() * A.transpose() * xyz_t;
		a1 = X0(0, 0);
		a2 = X0(1, 0);
		a3 = X0(2, 0);
		b1 = X0(3, 0);
		b2 = X0(4, 0);
		b3 = X0(5, 0);
		c1 = X0(6, 0);
		c2 = X0(7, 0);
		c3 = X0(8, 0);
		K = sqrt(1/3.0*(a1* a1 + a2 * a2 + a3 * a3+
			b1 * b1 + b2 * b2 + b3 * b3+
			c1 * c1 + c2 * c2 + c3 * c3));
		dx = X0(9, 0);
		dy = X0(10, 0);
		dz = X0(11, 0);
		R0 << a1, a2, a3,
			b1, b2, b3,
			c1, c2, c3;
		X << dx, dy, dz, K, a1, a2, a3, b1, b2, b3, c1, c2, c3;
		A.setZero(n, u);
	}
	//*********************//
	MatrixXd temp1(n, 1);//R*[x; y; z]
	MatrixXd temp2(n, 1);//[dx; dy; dz]
	MatrixXd B(s, u);//限制方程系数B s*u
	MatrixXd W(s, 1);//限制方程常数向量 s*1
	MatrixXd D(n + s, u);//设计矩阵(n+s)*u
	MatrixXd C(n + s, 1);//常数向量(n+s)*1
	
	x0.setZero(u, 1);
	x1.setOnes(u, 1);
	int count = 0;
	while((x1-x0).lpNorm<2>()>1e-10)
	{
		x1 = x0;
		count++;
		temp1.setZero(n, 1);
		temp1 << a1 * x_s + a2 * y_s + a3 * z_s,
			b1* x_s + b2 * y_s + b3 * z_s,
			c1* x_s + c2 * y_s + c3 * z_s;

		A << Kron(EYE3, ONEd), temp1, Kron(K * EYE3, xyz_s);
		temp2.setZero(n, 1);
		temp2 << dx * ONEd,
			dy* ONEd,
			dz* ONEd;
		l << temp2 + K * temp1 - xyz_t;
		B.setZero(s, 9);
		B << 2 * a1, 2 * a2, 2 * a3, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 2 * b1, 2 * b2, 2 * b3, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 2 * c1, 2 * c2, 2 * c3,
			a2, a1, 0, b2, b1, 0, c2, c1, 0,
			a3, 0, a1, b3, 0, b1, c3, 0, c1,
			0, a3, a2, 0, b3, b2, 0, c3, c2;
		MatrixXd tempB = B;
		B.setZero(s, u);
		B << MatrixXd::Zero(s, 4), tempB;
		
		W << a1 * a1 + a2 * a2 + a3 * a3 - 1,
			b1* b1 + b2 * b2 + b3 * b3 - 1,
			c1* c1 + c2 * c2 + c3 * c3 - 1,
			a1* a2 + b1 * b2 + c1 * c2,
			a1* a3 + b1 * b3 + c1 * c3,
			a2* a3 + b2 * b3 + c2 * c3;

		D << A,
			B;
		
		C << l,
			W;
		C = -C;
		//V=Dx-C
		//-e=Ax-y

		x0 = (D.transpose() * P * D).inverse() * D.transpose()* P * C;//求解法方程
		X = X + x0;
		dx = X(0, 0);
		dy = X(1, 0);
		dz = X(2, 0);
		K  = X(3, 0);
		a1 = X(4, 0); a2 = X(5, 0); a3 = X(6, 0);
		b1 = X(7, 0); b2 = X(8, 0); b3 = X(9, 0);
		c1 = X(10, 0); c2 = X(11, 0); c3 = X(12, 0);
		R0 << a1, a2, a3,
			b1, b2, b3,
			c1, c2, c3;
	}
	
	MatrixXd V(n, 1);//改正数
	
	V = (A * x0 - l);
	for (int i = 0; i < d; i++)
	{
		coorECEF_withname temp_e = {V(i,0), V(i + d,0), -V(i + 2 * d,0), name[i]};
		para.precision.e.push_back(temp_e);
	}
	//计算sigma0
	P.setIdentity(n, n);
	para.precision.sigma0 = ((V.transpose() *P* V) / para.precision.r)(0, 0);
	para.precision.sigma0 = sqrt(para.precision.sigma0);
	para.R = R0;
	para.K = K;
	para.dx = dx;
	para.dy = dy;
	para.dz = dz;
	para.count = count;
	return true;
}

bool computeDstCoorPara_13(vector<coorECEF_withname>src, vector<coorECEF_withname>& dst, parameter_13 para)
{
	int d = src.size();
	if (d == 0)
	{
		cout << "源坐标数据文件为空！\n";
		return false;
	}
	MatrixXd R(3, 3);
	R = para.R;
	double K = para.K;
	MatrixXd dxyz(3, 1);
	dxyz << para.dx,
		para.dy,
		para.dz;
	
	for (int i = 0; i < d; i++)
	{
		MatrixXd x_t(3, 1);
		MatrixXd x_s(3, 1);
		x_s << src[i].x,
			src[i].y,
			src[i].z;
		x_t = K * R * x_s + dxyz;
		coorECEF_withname temp = { x_t(0,0),x_t(1,0),x_t(2,0),src[i].name };
		dst.push_back(temp);
	}
	return true;
}

//计算残差
bool computeResidual(vector<coorECEF_withname>estimate, vector<coorECEF_withname> real, precisionResidual& res)
{
	map<string, coorECEF> src;
	map<string, coorECEF> dst;
	for (int i = 0; i < estimate.size(); i++)
	{
		coorECEF temp = { estimate[i].x,estimate[i].y,estimate[i].z };
		if (src.count(estimate[i].name) > 0)
		{
			cout << "检核点<估计值>数据内存在同名点！\n";
			return false;
		}
		else
			src.insert(pair<string,coorECEF>(estimate[i].name, temp));
	}
	for (int i = 0; i < real.size(); i++)
	{
		coorECEF temp = { real[i].x,real[i].y,real[i].z };
		if (dst.count(real[i].name) > 0)
		{
			cout << "检核点<真值>数据内存在同名点！\n";
			return false;
		}
		else
			dst.insert(pair<string, coorECEF>(real[i].name, temp));
	}
	for (auto it = src.begin(); it != src.end(); ++it)
	{
		if(dst.count(it->first)>0)
		{
			string name = it->first;
			coorECEF c1 = it->second;
			coorECEF c2 = dst[name];
			coorECEF_withname temp = { c2.x - c1.x, c2.y - c1.y, c2.z - c1.z, name };
			res.e.push_back(temp);
		}
		else
		{
			cout << "<估计值>数据与<真值>数据不匹配！\n<真值>数据中缺少该点: ";
			cout << it->first << " （点名）\n";
			return false;
		}
	}
	double variance = 0;
	for(int i=0;i<res.e.size();i++)
	{
		variance += (pow(res.e[i].x, 2) + pow(res.e[i].y, 2) + pow(res.e[i].x, 2));
	}
	res.variance = variance;
	return true;
}

//WGS84椭球参数
static double WGS84_a = 6378137;
static double WGS84_b = 6356752;
double WGS84_e2 = 0.00669437999013;

bool coorXYZ2BLH(vector<coorECEF_withname>src, vector<coorBLH_withname>&dst)
{
	int n = src.size();
	if (n == 0)
	{
		cout << "源坐标文件为空！\n";
		return false;
	}
	
	for (int i=0;i<n;i++)
	{
		double B0, B1, H;
		double x, y, z;
		x = src[i].x;
		y = src[i].y; 
		z = src[i].z;
		B0 = atan2(z , sqrt(x * x + y * y));
		B1 = B0 + 3;
		double N; 
		while(abs(B0-B1)>1e-6)
		{
			B1 = B0;
			B0 = atan(1/sqrt(x*x+y*y)*(z+WGS84_a*WGS84_e2*tan(B0)/sqrt(1+tan(B0)*tan(B0)-WGS84_e2*tan(B0)*tan(B0))));
		}
		coorBLH_withname temp;
		temp.name = src[i].name;
		temp.L=atan2(y, x);
		//东经正 西经负
		if (temp.L < 0 && y>0)
			temp.L = atan(1) * 4 + temp.L;
		else if(temp.L>0 && y<0)
			temp.L = temp.L - atan(1) * 4;

		N = WGS84_a / (sqrt(1 - WGS84_e2 * sin(B0) * sin(B0)));
		H = sqrt(x * x + y * y) / cos(B0) - N;
		temp.B = B0;
		temp.H = H;
		dst.push_back(temp);
	}
	return true;;
}

bool coorBLH2XYZ(vector<coorBLH_withname> src, vector<coorECEF_withname>& dst)
{
	int n = src.size();
	if (n == 0)
	{
		cout << "源坐标文件为空！\n";
		return false;
	}
	for (int i=0;i<n;i++)
	{
		double B, L, H;
		B = src[i].B;
		L = src[i].L;
		H = src[i].H;
		double N = WGS84_a / (sqrt(1 - WGS84_e2 * sin(B) * sin(B)));
		coorECEF_withname temp;
		temp.name = src[i].name;
		temp.x = (N + H) * cos(B) * cos(L);
		temp.y = (N + H) * cos(B) * sin(L);
		temp.z = (N * (1 - WGS84_e2) + H) * sin(B);
		dst.push_back(temp);
	}
}

bool coorXYZ2NEU(coorECEF_withname center0, vector<coorECEF_withname> src, vector<coorNEU_withname>& dst)
{
	vector<coorECEF_withname> temp1;
	vector<coorBLH_withname> temp2;
	temp1.push_back(center0);
	coorXYZ2BLH(temp1,temp2);
	MatrixXd mat(3, 3);
	coorBLH_withname center= temp2[0];
	
	mat(0, 0) = -sin(center.L);
	mat(0, 1) = cos(center.L);
	mat(0, 2) = 0;
	
	mat(1, 0) = -sin(center.B) * cos(center.L);
	mat(1, 1) = -sin(center.B) * sin(center.L);
	mat(1, 2) = cos(center.B);

	mat(2, 0) = cos(center.B) * cos(center.L);
	mat(2, 1) = cos(center.B) * sin(center.L);
	mat(2, 2) = sin(center.B);

	for(int i=0;i<src.size();i++)
	{
		double E, N, U;
		E = mat(0, 0) * (src[i].x - center0.x) + mat(0, 1) * (src[i].y - center0.y) + mat(0, 2) * (src[i].z - center0.z);
		N = mat(1, 0) * (src[i].x - center0.x) + mat(1, 1) * (src[i].y - center0.y) + mat(1, 2) * (src[i].z - center0.z);
		U = mat(2, 0) * (src[i].x - center0.x) + mat(2, 1) * (src[i].y - center0.y) + mat(2, 2) * (src[i].z - center0.z);
		coorNEU_withname temp = { N,E,U,src[i].name };
		dst.push_back(temp);
	}
	return true;
}