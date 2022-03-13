#pragma once
#include <fstream>
#include <vector>
#include <complex>
#include <iostream>

using std::vector;
using std::cout;
using std::complex;



bool compareDoubles(vector<double>& l, vector<double>& r)
	{
	if (l.size() != r.size()) return 0;
	if (l.empty() || r.empty()) return 0;
	auto it_r = r.begin();
	for (auto it_l = l.begin(); it_l != l.end(); it_l++, it_r++)
		{
		if (std::abs(*it_l - *it_r) > 0.1)
			{
			cout << "@@!= on " << *it_l << " - " << *it_r << std::endl;
			return 0;
			}

		}
	return 1;
	}