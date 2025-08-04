#pragma once

#include "Highs.h"
#include "ProblemData.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace bnp
{
	class RMP
	{
	public:
		Highs highs;
		HighsModel model;
		vector<int*> patterns;
		vector<int*> forbidden_patterns;

		RMP();
		~RMP();
		RMP(const RMP& other);               
		RMP& operator=(const RMP& other);  
		void copy_patterns_from(const RMP& other);
		void RMP::copy_forbidden_patterns_from(const RMP& other);

		void initialize();
		double solve_RMP(vector<double>& duals);
	};

	class SP
	{
	public:
		Highs highs;
		HighsModel SP_model;

		SP();
		bool solve_SP(RMP& rmp, const vector<double>& duals);
	};

	bool has_fractional_solution(RMP& rmp, int& branch_var, double& branch_value);

	bool is_forbidden_pattern(int* new_pat, const vector<int*>& forbidden_patterns, int n);
}