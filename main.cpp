// CoorTrans.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <sstream>
#include <ctime>
#include <fstream>
#include "CoorTrans.h"
using namespace std;
bool readECEFCoor(vector<coorECEF_withname>&);//读取ECEF坐标文件 文件格式：点名，X，Y，Z .xyz
bool readBLHCoor(vector<coorBLH_withname>& result);//读取BLH坐标文件 文件格式：点名，X，Y，Z .blh
bool writeBLHCoor(vector<coorBLH_withname> src, string &filename);//写入BLH坐标文件 文件格式：点名，B，L，H .txt
bool writeECEFCoor(vector<coorECEF_withname> src, string& filename);//写入BLH坐标文件 文件格式：点名，B，L，H .txt
bool writeNEUCoor(vector<coorNEU_withname> src, string& filename);//写入NEU坐标文件 文件格式：点名，N，E，U .txt
int getNumberPrecision(string x);//获取有效数字位数
int numberPrecision = 0;//有效数字位数
template <class Type>
Type stringToNum(const string& str)//将string转换为数字
{
	istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}
bool allIsInt(string str)//判断是否全是INT
{
	for (int i = 0; i < str.size(); i++)
	{
		int tmp = (int)str[i];
		if (tmp >= 48 && tmp <= 57)
		{
			continue;
		}
		else
		{
			return false;
		}
	}
	return true;
}
bool cmp(coorECEF_withname a, coorECEF_withname b)
{
	if(allIsInt(a.name)&&allIsInt(b.name))
	{
		return stringToNum<int>(a.name) < stringToNum<int>(b.name);
	}
	return a.name < b.name;
}
int main()
{
	cout << "----------------------------------" << endl;
	cout << "欢迎使用坐标转换程序！\n";
	cout << "https://github.com/AkexStar/CoorTrans/\n";
	cout << "Email: lijintao.alex@qq.com\n";
	cout << "----------------------------------" << endl;
	while (true)
	{
		cout << "\n--------请选择坐标转换类型！--------" << endl;
		cout << "<小角度四参数坐标转换>-------请输入1\n";
		cout << "<小角度六参数坐标转换>-------请输入2\n";
		cout << "<小角度七参数坐标转换>-------请输入3\n";
		cout << "<大角度十三参数坐标转换>-----请输入4\n";
		cout << "<地心地固系<->大地坐标系转换>请输入5\n";
		cout << "<地心地固系 ->站心坐标系转换>请输入6\n";
		cout << "<退出程序>-------------------请输入0\n";
		cout << "------------------------------------" << endl;
		
		string select_kind;
		cin >> select_kind;
		while(select_kind.size()!=1)
		{
			cout << "请输入0-6的提示数字选择下一步操作！\n" << endl;
			select_kind.clear();
			cin >> select_kind;
		}
		
		switch (stringToNum<int>(select_kind))
		{
		case 0:
			return 0;
		case 1://小角度四参数坐标转换
			{
				vector<coorECEF_withname> src_coor;	//源坐标
				vector<coorECEF_withname> dst_coor;	//目标坐标
				map<string, vector<coorECEF>> coor;	//坐标集合
				parameter_4 para;	//计算所得四参数
				getchar();	//读取多余的换行符
				
				cout << "读入源坐标系坐标：\n";
				if (readECEFCoor(src_coor))
				{
					for (int i = 0; i < src_coor.size(); i++)
					{
						coorECEF temp = { src_coor[i].x,src_coor[i].y,src_coor[i].z };
						if (coor.count(src_coor[i].name) > 0)
						{
							coor[src_coor[i].name].push_back(temp);
						}
						else
						{
							vector<coorECEF> temp2;
							temp2.push_back(temp);
							coor.insert(pair<string, vector<coorECEF>>(src_coor[i].name, temp2));
						}
					}
				}
				else break;
				
				cout << "读入目标坐标系坐标：\n";
				if (readECEFCoor(dst_coor))
				{
					for (int i = 0; i < dst_coor.size(); i++)
					{
						coorECEF temp = { dst_coor[i].x,dst_coor[i].y,dst_coor[i].z };
						if (coor.count(dst_coor[i].name) > 0)
						{
							coor[dst_coor[i].name].push_back(temp);
						}
						else
						{
							vector<coorECEF> temp2;
							temp2.push_back(temp);
							coor.insert(pair<string, vector<coorECEF>>(dst_coor[i].name, temp2));
						}
					}
				}
				else break;
				
				cout << "数据读入成功！下面进行计算参数！\n";
				if (computePara_4(coor, para))
				{
					cout << "-------------计算结果-------------\n";
					cout.setf(ios::showpoint); //将小数精度后面的0显示出来
					cout.precision(numberPrecision); //设置输出精度，保留有效数字
					cout << "条件方程数：\n\t" << para.precision.n << "\n" << "冗余度：\n\t" << para.precision.r << "\n";
					cout << "解得参数：\n";
					cout << "\tdx=" << para.dx;
					cout << "\n\tdy=" << para.dy;
					cout << "\n\tR=" << para.R;
					cout << "\n\tK=" << para.K;
					cout << "\n改正数：\n";
					for(int i=0;i<para.precision.e.size();i++)
					{
						//cout<<
						cout <<"\t"<< para.precision.e[i].name<<"\tVx="<<para.precision.e[i].x
						<<"\tVy="<<para.precision.e[i].y<<"\n";
					}
					cout << "Sigma0=\n\t"<<para.precision.sigma0 << endl;
					cout << "---------------------------------\n";
				}
				else break;

				//已知参数后进行坐标转换
				while(true)
				{
					cout << "是否继续在当前参数下求解目标坐标？\n是 请输入1 / 否 请输入0" << endl;
					int is_compute_int;
					cin >> is_compute_int;
					if (is_compute_int == 1)
					{
						cout << "读入需计算的源坐标系坐标：\n";
						vector<coorECEF_withname> src_compute_coor;//需要转换的源坐标
						vector<coorECEF_withname> dst_compute_coor;//转换所得的目标坐标
						getchar();
						if (readECEFCoor(src_compute_coor))
						{
							computeDstCoorPara_4(src_compute_coor, dst_compute_coor, para);
						}
						else
							break;
						//cout.setf(ios::showpoint); //将小数精度后面的0显示出来
						cout.flags(ios::fixed);
						cout.precision(3); //设置输出精度，保留有效数字
						cout << "-------------计算结果-------------\n";
						cout << "已转换完" << dst_compute_coor.size() << "个坐标：\n";
						cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
						for(int i=0;i<dst_compute_coor.size();i++)
						{
							cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].x
								<< "\t" << dst_compute_coor[i].y << "\t" << dst_compute_coor[i].z << endl;
						}
						cout << "---------------------------------\n";
						cout << "是否对当前已计算得到的坐标<估计值>进行外检核？\n是 请输入1 / 否 请输入0" << endl;
						int is_judge_int;
						cin >> is_judge_int;
						if (is_judge_int == 1)
						{
							cout << "读入检核点<真值>坐标：\n";
							getchar();
							vector<coorECEF_withname> real_judge_coor;//真值坐标
							precisionResidual res;
							if (readECEFCoor(real_judge_coor))
							{
								computeResidual(dst_compute_coor, real_judge_coor, res);
							}
							else
								break;
							cout << "-------------计算结果-------------\n";
							cout.flags(ios::fixed);
							cout.precision(4); //设置保留小数位数
							cout << "残差平方和：\n\t" << res.variance << "\n残差：\n";
							sort(res.e.begin(), res.e.end(), cmp);
							cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
							for (int i = 0; i < res.e.size(); i++)
							{
								cout << "\t" << res.e[i].name << "\t" << res.e[i].x
									<< "\t" << res.e[i].y << "\t" << res.e[i].z << endl;
							}
							cout << "---------------------------------\n";
						}
					}
					else if(is_compute_int==0)
						break;
				}
				break;
			}
		case 2://小角度六参数坐标转换
			{
				vector<coorECEF_withname> src_coor;	//源坐标
				vector<coorECEF_withname> dst_coor;	//目标坐标
				map<string, vector<coorECEF>> coor;	//坐标集合
				parameter_6 para;	//计算所得七参数
				getchar();	//读取多余的换行符

				cout << "读入源坐标系坐标：\n";
				if (readECEFCoor(src_coor))
				{
					for (int i = 0; i < src_coor.size(); i++)
					{
						coorECEF temp = { src_coor[i].x,src_coor[i].y,src_coor[i].z };
						if (coor.count(src_coor[i].name) > 0)
						{
							coor[src_coor[i].name].push_back(temp);
						}
						else
						{
							vector<coorECEF> temp2;
							temp2.push_back(temp);
							coor.insert(pair<string, vector<coorECEF>>(src_coor[i].name, temp2));
						}
					}
				}
				else break;

				cout << "读入目标坐标系坐标：\n";
				if (readECEFCoor(dst_coor))
				{
					for (int i = 0; i < dst_coor.size(); i++)
					{
						coorECEF temp = { dst_coor[i].x,dst_coor[i].y,dst_coor[i].z };
						if (coor.count(dst_coor[i].name) > 0)
						{
							coor[dst_coor[i].name].push_back(temp);
						}
						else
						{
							vector<coorECEF> temp2;
							temp2.push_back(temp);
							coor.insert(pair<string, vector<coorECEF>>(dst_coor[i].name, temp2));
						}
					}
				}
				else break;

				cout << "数据读入成功！下面进行计算参数！\n";
				if (computePara_6(coor, para))
				{
					cout << "-------------计算结果-------------\n";
					cout.setf(ios::showpoint); //将小数精度后面的0显示出来
					cout.precision(numberPrecision); //设置输出精度，保留有效数字
					cout << "条件方程数：\n\t" << para.precision.n << "\n" << "冗余度：\n\t" << para.precision.r << "\n";
					cout << "解得参数：\n";
					cout << "\tdx=" << para.dx;
					cout << "\n\tdy=" << para.dy;
					cout << "\n\tdz=" << para.dz;
					cout << "\n\tRx=" << para.Rx;
					cout << "\n\tRy=" << para.Ry;
					cout << "\n\tRz=" << para.Rz;
					cout << "\n改正数：\n";
					for (int i = 0; i < para.precision.e.size(); i++)
					{
						//cout<<
						cout << "\t" << para.precision.e[i].name << "\tVx=" << para.precision.e[i].x
							<< "\tVy=" << para.precision.e[i].y << "\tVz=" << para.precision.e[i].z << "\n";
					}
					cout << "Sigma0=\n\t" << para.precision.sigma0 << endl;
					cout << "---------------------------------\n";
				}
				else break;

				//已知参数后进行坐标转换
				while (true)
				{
					cout << "是否继续在当前参数下求解目标坐标？\n是 请输入1 / 否 请输入0" << endl;
					int is_compute_int;
					cin >> is_compute_int;
					if (is_compute_int == 1)
					{
						cout << "读入需计算的源坐标系坐标：\n";
						vector<coorECEF_withname> src_compute_coor;//需要转换的源坐标
						vector<coorECEF_withname> dst_compute_coor;//转换所得的目标坐标
						getchar();
						if (readECEFCoor(src_compute_coor))
						{
							computeDstCoorPara_6(src_compute_coor, dst_compute_coor, para);
						}
						else
							break;
						//cout.setf(ios::showpoint); //将小数精度后面的0显示出来
						cout.flags(ios::fixed);
						cout.precision(4); //设置保留小数位数
						cout << "-------------计算结果-------------\n";
						cout << "已转换完" << dst_compute_coor.size() << "个坐标：\n";
						cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
						for (int i = 0; i < dst_compute_coor.size(); i++)
						{
							cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].x
								<< "\t" << dst_compute_coor[i].y << "\t" << dst_compute_coor[i].z << endl;
						}
						cout << "---------------------------------\n";
						cout << "是否对当前已计算得到的坐标<估计值>进行外检核？\n是 请输入1 / 否 请输入0" << endl;
						int is_judge_int;
						cin >> is_judge_int;
						if (is_judge_int == 1)
						{
							cout << "读入检核点<真值>坐标：\n";
							getchar();
							vector<coorECEF_withname> real_judge_coor;//真值坐标
							precisionResidual res;
							if (readECEFCoor(real_judge_coor))
							{
								computeResidual(dst_compute_coor, real_judge_coor, res);
							}
							else
								break;
							cout << "-------------计算结果-------------\n";
							cout.flags(ios::fixed);
							cout.precision(4); //设置保留小数位数
							cout << "残差平方和：\n\t" << res.variance << "\n残差：\n";
							sort(res.e.begin(), res.e.end(), cmp);
							cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
							for (int i = 0; i < res.e.size(); i++)
							{
								cout << "\t" << res.e[i].name << "\t" << res.e[i].x
									<< "\t" << res.e[i].y << "\t" << res.e[i].z << endl;
							}
							cout << "---------------------------------\n";
						}
					}
					else if (is_compute_int == 0)
						break;
				}
				break;
			}
		case 3://小角度七参数坐标转换
			{
				vector<coorECEF_withname> src_coor;	//源坐标
				vector<coorECEF_withname> dst_coor;	//目标坐标
				map<string, vector<coorECEF>> coor;	//坐标集合
				parameter_7 para;	//计算所得七参数
				getchar();	//读取多余的换行符

				cout << "读入源坐标系坐标：\n";
				if (readECEFCoor(src_coor))
				{
					for (int i = 0; i < src_coor.size(); i++)
					{
						coorECEF temp = { src_coor[i].x,src_coor[i].y,src_coor[i].z };
						if (coor.count(src_coor[i].name) > 0)
						{
							coor[src_coor[i].name].push_back(temp);
						}
						else
						{
							vector<coorECEF> temp2;
							temp2.push_back(temp);
							coor.insert(pair<string, vector<coorECEF>>(src_coor[i].name, temp2));
						}
					}
				}
				else break;

				cout << "读入目标坐标系坐标：\n";
				if (readECEFCoor(dst_coor))
				{
					for (int i = 0; i < dst_coor.size(); i++)
					{
						coorECEF temp = { dst_coor[i].x,dst_coor[i].y,dst_coor[i].z };
						if (coor.count(dst_coor[i].name) > 0)
						{
							coor[dst_coor[i].name].push_back(temp);
						}
						else
						{
							vector<coorECEF> temp2;
							temp2.push_back(temp);
							coor.insert(pair<string, vector<coorECEF>>(dst_coor[i].name, temp2));
						}
					}
				}
				else break;

				cout << "数据读入成功！下面进行计算参数！\n";
				if (computePara_7(coor, para))
				{
					cout << "-------------计算结果-------------\n";
					cout.setf(ios::showpoint); //将小数精度后面的0显示出来
					cout.precision(numberPrecision); //设置输出精度，保留有效数字
					cout << "条件方程数：\n\t" << para.precision.n << "\n" << "冗余度：\n\t" << para.precision.r << "\n";
					cout << "解得参数：\n";
					cout << "\tdx=" << para.dx;
					cout << "\n\tdy=" << para.dy;
					cout << "\n\tdz=" << para.dz;
					cout << "\n\tRx=" << para.Rx;
					cout << "\n\tRy=" << para.Ry;
					cout << "\n\tRz=" << para.Rz;
					cout << "\n\tK=" << para.K;
					cout << "\n改正数：\n";
					for (int i = 0; i < para.precision.e.size(); i++)
					{
						//cout<<
						cout << "\t" << para.precision.e[i].name << "\tVx=" << para.precision.e[i].x
							<< "\tVy=" << para.precision.e[i].y << "\tVz=" << para.precision.e[i].z << "\n";
					}
					cout << "Sigma0=\n\t" << para.precision.sigma0 << endl;
					cout << "---------------------------------\n";
				}
				else break;

				//已知参数后进行坐标转换
				while (true)
				{
					cout << "是否继续在当前参数下求解目标坐标？\n是 请输入1 / 否 请输入0" << endl;
					int is_compute_int;
					cin >> is_compute_int;
					if (is_compute_int == 1)
					{
						cout << "读入需计算的源坐标系坐标：\n";
						vector<coorECEF_withname> src_compute_coor;//需要转换的源坐标
						vector<coorECEF_withname> dst_compute_coor;//转换所得的目标坐标
						getchar();
						if (readECEFCoor(src_compute_coor))
						{
							computeDstCoorPara_7(src_compute_coor, dst_compute_coor, para);
						}
						else
							break;
						//cout.setf(ios::showpoint); //将小数精度后面的0显示出来
						cout.flags(ios::fixed);
						cout.precision(4); //设置保留小数位数
						cout << "-------------计算结果-------------\n";
						cout << "已转换完" << dst_compute_coor.size() << "个坐标：\n";
						cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
						for (int i = 0; i < dst_compute_coor.size(); i++)
						{
							cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].x
								<< "\t" << dst_compute_coor[i].y << "\t" << dst_compute_coor[i].z << endl;
						}
						cout << "---------------------------------\n";
						{
							cout << "是否对当前已计算得到的坐标<估计值>进行外检核？\n是 请输入1 / 否 请输入0" << endl;
							int is_judge_int;
							cin >> is_judge_int;
							if (is_judge_int == 1)
							{
								cout << "读入检核点<真值>坐标：\n";
								getchar();
								vector<coorECEF_withname> real_judge_coor;//真值坐标
								precisionResidual res;
								if (readECEFCoor(real_judge_coor))
								{
									computeResidual(dst_compute_coor, real_judge_coor, res);
								}
								else
									break;
								cout << "-------------计算结果-------------\n";
								cout.flags(ios::fixed);
								cout.precision(4); //设置保留小数位数
								cout << "残差平方和：\n\t" << res.variance << "\n残差：\n";
								sort(res.e.begin(), res.e.end(), cmp);
								cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
								for (int i = 0; i < res.e.size(); i++)
								{
									cout << "\t" << res.e[i].name << "\t" << res.e[i].x
										<< "\t" << res.e[i].y << "\t" << res.e[i].z << endl;
								}
								cout << "---------------------------------\n";
							}
						}
					}
					else if (is_compute_int == 0)
						break;
				}
				break;
			}
		case 4://大角度十三参数坐标转换
			{
			vector<coorECEF_withname> src_coor;	//源坐标
			vector<coorECEF_withname> dst_coor;	//目标坐标
			map<string, vector<coorECEF>> coor;	//坐标集合
			parameter_13 para;	//计算所得十三参数
			getchar();	//读取多余的换行符

			cout << "读入源坐标系坐标：\n";
			if (readECEFCoor(src_coor))
			{
				for (int i = 0; i < src_coor.size(); i++)
				{
					coorECEF temp = { src_coor[i].x,src_coor[i].y,src_coor[i].z };
					if (coor.count(src_coor[i].name) > 0)
					{
						coor[src_coor[i].name].push_back(temp);
					}
					else
					{
						vector<coorECEF> temp2;
						temp2.push_back(temp);
						coor.insert(pair<string, vector<coorECEF>>(src_coor[i].name, temp2));
					}
				}
			}
			else break;

			cout << "读入目标坐标系坐标：\n";
			if (readECEFCoor(dst_coor))
			{
				for (int i = 0; i < dst_coor.size(); i++)
				{
					coorECEF temp = { dst_coor[i].x,dst_coor[i].y,dst_coor[i].z };
					if (coor.count(dst_coor[i].name) > 0)
					{
						coor[dst_coor[i].name].push_back(temp);
					}
					else
					{
						vector<coorECEF> temp2;
						temp2.push_back(temp);
						coor.insert(pair<string, vector<coorECEF>>(dst_coor[i].name, temp2));
					}
				}
			}
			else break;

			cout << "数据读入成功！下面进行计算参数！\n";
			if (computePara_13(coor, para))
			{
				cout << "-------------计算结果-------------\n";
				cout.setf(ios::showpoint); //将小数精度后面的0显示出来
				cout.precision(numberPrecision); //设置输出精度，保留有效数字
				cout << "条件方程数：\n\t" << para.precision.n << "\n" << "冗余度：\n\t" << para.precision.r << "\n";
				cout << "解得参数：\n";
				cout << "\tdx=" << para.dx;
				cout << "\n\tdy=" << para.dy;
				cout << "\n\tdz=" << para.dz;
				cout << "\n\tR=\n" << para.R;
				cout << "\n\tK=" << para.K;
				cout << "\n迭代次数=\n\t" << para.count;
				cout << "\n改正数：\n";
				cout << "\tPointName\tVx\tVy\tVz\n";
				sort(para.precision.e.begin(), para.precision.e.end(), cmp);
				for (int i = 0; i < para.precision.e.size(); i++)
				{
					cout << "\t" << para.precision.e[i].name << "\t" << para.precision.e[i].x
						<< "\t" << para.precision.e[i].y << "\t" << para.precision.e[i].z << "\n";
				}
				cout << "Sigma0=\n\t" << para.precision.sigma0 << endl;
				cout << "---------------------------------\n";
				
			}
			else break;

			//已知参数后进行坐标转换
			while (true)
			{
				cout << "是否继续在当前参数下求解目标坐标？\n是 请输入1 / 否 请输入0" << endl;
				int is_compute_int;
				cin >> is_compute_int;
				if (is_compute_int == 1)
				{
					cout << "读入需计算的源坐标系坐标：\n";
					vector<coorECEF_withname> src_compute_coor;//需要转换的源坐标
					vector<coorECEF_withname> dst_compute_coor;//转换所得的目标坐标
					getchar();
					if (readECEFCoor(src_compute_coor))
					{
						computeDstCoorPara_13(src_compute_coor, dst_compute_coor, para);
					}
					else
						break;
					//cout.setf(ios::showpoint); //将小数精度后面的0显示出来
					cout.flags(ios::fixed);
					cout.precision(numberPrecision); //设置保留小数位数
					cout << "-------------计算结果-------------\n";
					cout << "已转换完" << dst_compute_coor.size() << "个坐标：\n";
					cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
					for (int i = 0; i < dst_compute_coor.size(); i++)
					{
						cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].x
							<< "\t" << dst_compute_coor[i].y << "\t" << dst_compute_coor[i].z << endl;
					}
					cout << "---------------------------------\n";
					cout << "是否对当前已计算得到的坐标<估计值>进行外检核？\n是 请输入1 / 否 请输入0" << endl;
					int is_judge_int;
					cin >> is_judge_int;
					if(is_judge_int==1)
					{
						cout << "读入检核点<真值>坐标：\n";
						getchar();
						vector<coorECEF_withname> real_judge_coor;//真值坐标
						precisionResidual res;
						if (readECEFCoor(real_judge_coor))
						{
							computeResidual(dst_compute_coor, real_judge_coor, res);
						}
						else
							break;
						cout << "-------------计算结果-------------\n";
						cout.flags(ios::fixed);
						cout.precision(4); //设置保留小数位数
						cout << "残差平方和：\n\t"<<res.variance<<"\n残差：\n";
						sort(res.e.begin(), res.e.end(),cmp);
						cout << "\tPointName\t" << "X\t" << "Y\t" << "Z" << endl;
						for(int i=0;i<res.e.size();i++)
						{
							cout << "\t" << res.e[i].name << "\t" << res.e[i].x
								<< "\t" << res.e[i].y << "\t" << res.e[i].z << endl;
						}
						cout << "---------------------------------\n";
					}
				}
				else if (is_compute_int == 0)
					break;
			}
				break;
			}
		case 5://地心地固系-大地坐标系转换
			{
				getchar();	//读取多余的换行符
				cout << "------------------------------------" << endl;
				cout << "ECEF_xyz -> BLH 请输入 1\n" << "BLH -> ECEF_xyz 请输入 2\n";
				int transform_selected;
				cin >> transform_selected;
				if(transform_selected==1)
				{
					cout << "读入源坐标系坐标：\n";
					vector<coorECEF_withname> src_coor;	//源坐标
					vector<coorBLH_withname> dst_compute_coor;	//目标坐标
					getchar();
					if (readECEFCoor(src_coor))
					{
						cout << "数据读入成功！下面进行计算！\n";
						if (coorXYZ2BLH(src_coor, dst_compute_coor))
						{
							/*cout << "-------------计算结果-------------\n";
							cout << "已转换完" << dst_compute_coor.size() << "个坐标：\n";
							cout << "\tPointName\t" << "B\t" << "L\t" << "H" << endl;
							for(int i=0;i<dst_compute_coor.size();i++)
							{
								cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].B
									<< "\t" << dst_compute_coor[i].L << "\t" << dst_compute_coor[i].H << endl;
							}
							cout << "---------------------------------\n";*/
							string filename="";
							if(writeBLHCoor(dst_compute_coor,filename))
							{
								cout << "数据文件已经保存至\n";
								cout << filename << endl;
							}
							cout << "---------------------------------\n";
						}
					}
					else break;
				}
				else if(transform_selected==2)
				{
					cout << "读入源坐标系坐标：\n";
					vector<coorBLH_withname> src_coor;	//源坐标
					vector<coorECEF_withname> dst_compute_coor;	//目标坐标
					getchar();
					if (readBLHCoor(src_coor))
					{
						cout << "数据读入成功！下面进行计算！\n";
						if (coorBLH2XYZ(src_coor, dst_compute_coor))
						{
							/*cout << "-------------计算结果-------------\n";
							cout << "已转换完" << dst_compute_coor.size() << "个坐标：\n";
							cout << "\tPointName\t" << "B\t" << "L\t" << "H" << endl;
							for(int i=0;i<dst_compute_coor.size();i++)
							{
								cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].B
									<< "\t" << dst_compute_coor[i].L << "\t" << dst_compute_coor[i].H << endl;
								cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].x
									<< "\t" << dst_compute_coor[i].y << "\t" << dst_compute_coor[i].z << endl;
							}
							cout << "---------------------------------\n";*/
							string filename = "";
							if (writeECEFCoor(dst_compute_coor, filename))
							{
								cout << "数据文件已经保存至\n";
								cout << filename << endl;
							}
							cout << "------------------------------------\n";
						}
					}
					else break;
				}
				break;
			}
		case 6://地心地固系-站心坐标系转换
			{
				
				cout << "------------------------------------" << endl;
				
				cout << "读入源坐标系坐标：\n";
				vector<coorECEF_withname> src_coor;	//源坐标
				vector<coorNEU_withname> dst_compute_coor;	//目标坐标
				getchar();	//读取多余的换行符
				if (readECEFCoor(src_coor))
				{
					cout << "数据读入成功！下面进行计算！\n";
					coorECEF_withname center = { src_coor[0].x,src_coor[0].y,src_coor[0].z,src_coor[0].name};
					if (coorXYZ2NEU(center,src_coor, dst_compute_coor))
					{
						/*cout << "-------------计算结果-------------\n";
						cout << "已转换完" << dst_compute_coor.size() << "个坐标：\n";
						cout << "\tPointName\t" << "B\t" << "L\t" << "H" << endl;
						for(int i=0;i<dst_compute_coor.size();i++)
						{
							cout << "\t" << dst_compute_coor[i].name << "\t" << dst_compute_coor[i].B
								<< "\t" << dst_compute_coor[i].L << "\t" << dst_compute_coor[i].H << endl;
						}
						cout << "---------------------------------\n";*/
						string filename = "";
						if (writeNEUCoor(dst_compute_coor, filename))
						{
							cout << "数据文件已经保存至\n";
							cout << filename << endl;
						}
						cout << "------------------------------------\n";
					}
				}
				break;
			}
		default:
			{
				cout << ">>输入有误！请输入0-6的数字进行菜单选择！\n<<";
				break;
			}
		}
	}
	return 0;
}

bool readECEFCoor(vector<coorECEF_withname> &result)
{
	string fileName;//文件名
	cout << "请输入坐标文件的路径！\n";
	getline(cin,fileName);
	int length = fileName.length();
	for(int i=0;i<length;i++)
	{
		if (fileName[i] == '\\')
		{
			fileName.insert(i, "\\");
			i++;
		}
	}
	fstream file;
	file.open(fileName);//读入文件
	if (!file)
	{
		cout << "文件地址错误！打开失败！" << endl;
		return false;
	}
	string coorDateStr;	//临时储存数据文件的每一行
	numberPrecision = 0;
	//一行一行读入
	while (getline(file, coorDateStr))
	{
		for (int i = 0; i < coorDateStr.length(); i++)
		{
			if (coorDateStr[i] == ',')
				coorDateStr[i] = '\t';
		}

		istringstream temp(coorDateStr);
		string Name;
		double X = 0, Y = 0, Z = 0;
		string x_string;
		string y_string;
		string z_string;
		temp >> Name;
		temp >> x_string;
		temp >> y_string;
		temp >> z_string;
		X = stringToNum<double>(x_string);
		Y = stringToNum<double>(y_string);
		Z = stringToNum<double>(z_string);
		numberPrecision = max(numberPrecision,(getNumberPrecision(x_string),
			max(getNumberPrecision(y_string), getNumberPrecision(z_string))));
		//判断有效数字位数
		coorECEF_withname tempCoor = {X, Y, Z, Name};
		result.push_back(tempCoor);
	}
	if (result.empty())
	{
		cout << "数据文件为空！请检查文件！" << endl;
		return false;
	}
		
	file.close();
	return true;
}

bool readBLHCoor(vector<coorBLH_withname>& result)
{
	string fileName;//文件名
	cout << "请输入坐标文件的路径！\n";
	getline(cin, fileName);
	int length = fileName.length();
	for (int i = 0; i < length; i++)
	{
		if (fileName[i] == '\\')
		{
			fileName.insert(i, "\\");
			i++;
		}
	}
	fstream file;
	file.open(fileName);//读入文件
	if (!file)
	{
		cout << "文件地址错误！打开失败！" << endl;
		return false;
	}
	string coorDateStr;	//临时储存数据文件的每一行
	numberPrecision = 0;
	//一行一行读入
	while (getline(file, coorDateStr))
	{
		for (int i = 0; i < coorDateStr.length(); i++)
		{
			if (coorDateStr[i] == ',')
				coorDateStr[i] = '\t';
		}

		istringstream temp(coorDateStr);
		string Name;
		double X = 0, Y = 0, Z = 0;
		string x_string;
		string y_string;
		string z_string;
		temp >> Name;
		temp >> x_string;
		temp >> y_string;
		temp >> z_string;
		X = stringToNum<double>(x_string);
		Y = stringToNum<double>(y_string);
		Z = stringToNum<double>(z_string);
		numberPrecision = max(numberPrecision, (getNumberPrecision(x_string),
			max(getNumberPrecision(y_string), getNumberPrecision(z_string))));
		//判断有效数字位数
		coorBLH_withname tempCoor = { X, Y, Z, Name };
		result.push_back(tempCoor);
	}
	if (result.empty())
	{
		cout << "数据文件为空！请检查文件！" << endl;
		return false;
	}

	file.close();
	return true;
}

bool writeBLHCoor(vector<coorBLH_withname> src, string &filename)
{
	if(filename=="")
	{
		time_t now = time(0);
		filename = ctime(&now);
		filename.pop_back();
		for (int i = 0; i < filename.size(); i++)
			if (filename[i] == ' '|| filename[i] == ':')
				filename[i] = '_';
		filename = "D:\\result_BLH_"+filename +".txt";
	}
	ofstream file;
	file.open(filename);
	if(file)
	{
		file.flags(ios::fixed);
		for(int i=0;i<src.size();i++)
		{
			src[i].B = src[i].B*180.0 / (atan(1) * 4);
			src[i].L = src[i].L*180.0 / (atan(1) * 4);
			//转换为度
			cout.flags(ios::fixed);
			file.precision(9); //设置保留小数位数
			file << src[i].name << "," << src[i].B << "," << src[i].L << ",";
			file.precision(3); //设置保留小数位数
			file << src[i].H << endl;
		}
		file.close();
		return true;
	}
	return false;
	
}

bool writeECEFCoor(vector<coorECEF_withname> src, string& filename)
{
	if (filename == "")
	{
		time_t now = time(0);
		filename = ctime(&now);
		filename.pop_back();
		for (int i = 0; i < filename.size(); i++)
			if (filename[i] == ' ' || filename[i] == ':')
				filename[i] = '_';
		filename = "D:\\result_BLH_" + filename + ".txt";
	}
	ofstream file;
	file.open(filename);
	if (file)
	{
		file.flags(ios::fixed);
		for (int i = 0; i < src.size(); i++)
		{
			cout.flags(ios::fixed);
			file.precision(9); //设置保留小数位数
			file << src[i].name << "," << src[i].x << "," << src[i].y << ","<< src[i].z << endl;
		}
		file.close();
		return true;
	}
	return false;

}

bool writeNEUCoor(vector<coorNEU_withname> src, string& filename)
{
	if (filename == "")
	{
		time_t now = time(0);
		filename = ctime(&now);
		filename.pop_back();
		for (int i = 0; i < filename.size(); i++)
			if (filename[i] == ' ' || filename[i] == ':')
				filename[i] = '_';
		filename = "D:\\result_NEU_" + filename + ".txt";
	}
	ofstream file;
	file.open(filename);
	if (file)
	{
		file.flags(ios::fixed);
		for (int i = 0; i < src.size(); i++)
		{
			cout.flags(ios::fixed);
			file.precision(9); //设置保留小数位数
			file << src[i].name << "," << src[i].N << "," << src[i].E << "," << src[i].U << endl;
		}
		file.close();
		return true;
	}
	return false;
}

int getNumberPrecision(string x)
{
	if (x.find('.') != string::npos)
		return x.size() - 1;
	return x.size();
}
