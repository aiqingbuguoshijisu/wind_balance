﻿#pragma once
#include<fstream>
#include<vector>
#include<string>
#include<regex>
#include<libxl.h>
#include<iostream>
#include<sstream>
#include<Eigen/Dense>
#include<Eigen/Core>

using namespace std;
using namespace Eigen;

wstring s2ws(const string& str)//将标准字符串（string）转换为宽字符字符串（wstring）这个是用来解决路径数据类型问题
{
	return wstring(str.begin(), str.end());
}

vector<vector<int>> find_zeroRows_index(const string&filePath)
{
	/*找到loadxx文件中载荷为0的行的索引
	* 前72行中，每12行中的载荷为0的行的索引放到一个vector中
	* 后72行中，每18行中的载荷为0的行的索引放到一个vector中
	*/
	vector<vector<int>> zero_rows;//存放载荷为0的行的索引

	ifstream inFile(filePath);//打开文件
	if (!inFile.is_open())
	{
		cerr << "Unable to open file for reading: " << filePath << endl;
		return zero_rows;
	}

	vector<string> lines;//将文件中内容读入lines数组，后面修改内容，也是直接会在这个数组中修改
	string line;//临时存放文件的一行内容，属于临时变量
	while (getline(inFile, line))
	{
		lines.push_back(line);
	}
	inFile.close();

	// 删除前三行
	lines.erase(lines.begin(), lines.begin() + 3);//loadxx文件中前三行没有

	regex pattern(R"(^\s*\d{1,3}\s+)"); //匹配每一行的行号
	for (auto& l : lines)
	{
		l = regex_replace(l, pattern, "");//删除每一行的行号
	}

	auto check_zero_row = [](const string& line) -> bool 
		{
			/*
			检查某一行的内容中是否全为0
			*/
		istringstream iss(line);//将字符串流line读入iss
		double value;
		while (iss >> value) //这个循环会将字符串流line中的每一个数值读入value，判断是否为0
		{
			if (value != 0) 
			{
				return false;
			}
		}
		return true;//如果全是0，则返回true，否则返回false
		};

	for (int i = 0; i < 72; i += 12)//处理前72行的载荷为0的行
	{
		vector<int> zero_row;
		for (int j = 0; j < 12; j++)
		{
			int index = i + j;
			if (check_zero_row(lines[index]))
			{
				zero_row.push_back(index);
			}
		}
		zero_rows.push_back(zero_row);
	}

	for (int i = 72; i < 144; i += 18)//处理后72行的载荷为0的行
	{
		vector<int> zero_row;
		for (int j = 0; j < 18; j++)
		{
			int index = i + j;
			if (check_zero_row(lines[index]))
			{
				zero_row.push_back(index);
			}
		}
		zero_rows.push_back(zero_row);
	}

	return zero_rows;
}

MatrixXd compute_cofficient_Y(const string& dataxxFilePath,const string& loadxxFilePath)
{
	ifstream inFile(dataxxFilePath);

	if (!inFile.is_open())
	{
		cerr << "Unable to open file for reading: " << dataxxFilePath << endl;
		return MatrixXd();
	}

	vector<string> lines;
	string line;

	while (getline(inFile, line))//将文件的每一line读入lines数组
	{
		lines.push_back(line);
	}
	inFile.close();

	if (lines.empty())
	{
		cerr << "File is empty or failed to read lines: " << dataxxFilePath << endl;
		return MatrixXd();
	}

	lines.erase(lines.begin());//删除lines数组的第一行内容

	regex pattern1(R"(^\s*\d{1,3}\s+)"); //匹配每一行的行号
	for (auto& l : lines)
	{
		l = regex_replace(l, pattern1, "");//删除每一行的行号
	}

	regex pattern2(R"((-?\d+\.\d+e[+-]\d+|\d+\.\d+|\d+)\s+)");//提取dataxx文件中的数据
	smatch matches;

	int rowCount = lines.size();//144
	int colCount = 12;//只要dataxx的前12列，其中前6列为电压值，后6列为公斤值，需要转换成牛
	MatrixXd dataxxx(lines.size(), colCount);//存放数据

	int rowIndex = 0;

	for (const auto& l : lines)
	{
		/*
		* 这个循环就是在dataxx文件中每提取一个数据，就直接放到dataxxx矩阵中
		*/
		int colIndex = 0;
		string::const_iterator searchStart(l.cbegin());//字符串迭代器，指向当前搜索的起始位置，每次迭代后会指向下一个匹配位置
		while (regex_search(searchStart, l.cend(), matches, pattern2) && colIndex < colCount)//直接使用line.cend()作为结束位置，是因为它不用动
		{
			dataxxx(rowIndex, colIndex) = stod(matches[0]);//将提取的字符串转换为double类型，string to double
			searchStart = matches.suffix().first;
			colIndex++;
		}
		if (colIndex == colCount) 
		{
			rowIndex++;  // 只有当这一行的数据完整时，才计入行数
		}
		else 
		{
			cerr << "Row does not contain exactly 12 values, skipping: " << l << endl;
		}
	}

	vector<MatrixXd> f_list;
	MatrixXd f = dataxxx.block(0, 6, rowCount, 6) * 9.8035;
	for (int i = 0; i < 6; i++)
	{
		f_list.push_back(f.col(i));
	}

	MatrixXd u = dataxxx.block(0, 0, rowCount, 6);

	MatrixXd compute_u_0 = u.col(1) - u.col(0);

	MatrixXd u_0 = MatrixXd::Zero(rowCount, 1);

	vector<vector<int>> zero_rows = find_zeroRows_index(loadxxFilePath);

	vector<double> zero_rows_u_0(zero_rows.size(),0.0);
	for (int i = 0; i < zero_rows.size(); i++)
	{
		double sum = 0;
		for (int j : zero_rows[i])
		{
			sum += compute_u_0(j,0);
		}
		zero_rows_u_0[i] = sum / zero_rows[i].size();
	}

	for (int k = 0; k < 72; k += 12) 
	{
		int i = k / 12;
		u_0.block(12 * i, 0, 12, 1) = compute_u_0.block(12 * i, 0, 12, 1).array() - zero_rows_u_0[i];
	}

	for (int k = 72; k < 144; k += 18) 
	{
		int i = (k - 72) / 18;
		u_0.block(72 + 18 * i, 0, 18, 1) = compute_u_0.block(72 + 18 * i, 0, 18, 1).array() - zero_rows_u_0[i + 6];
	}

	MatrixXd py_one_infer(rowCount, 5);
	py_one_infer << f.col(1), f.col(2), f.col(3), f.col(4), f.col(5);

	MatrixXd py_two_infer(rowCount, 21);
	for (int i = 0; i < 6; i++)
	{
		py_two_infer.col(i) = f.col(i).array().square();
	}

	int index = 6;
	for (int i = 0; i < 5; i++)
	{
		for (int j = i + 1; j < 6; j++)
		{
			py_two_infer.col(index++) = f.col(i).array() * f.col(j).array();
		}
	}

	MatrixXd A(144, 27);
	A<<u_0, py_one_infer, py_two_infer;

	MatrixXd X = (A.transpose() * A).inverse() * A.transpose() * f.col(0);

	return X;
}