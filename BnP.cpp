#include "BnP.h"
#include "Highs.h"
#include "ProblemData.h"
#include <limits>
#include <iostream>
#include <cmath>

using namespace std;

namespace bnp
{
	RMP::RMP()
	{
		model.lp_.sense_ = ObjSense::kMinimize;
		highs.setOptionValue("output_flag", false);
	}

	RMP::~RMP()
	{
		for (auto p : patterns)
		{
			delete[] p;
		}

		for (auto p : forbidden_patterns)
		{
			delete[] p;
		}
	}

	// 복사 생성자
	RMP::RMP(const RMP& other)
	{
		model = other.model;
		highs.passModel(model);  // HiGHS는 반드시 모델 전달 필요
		copy_patterns_from(other);
		copy_forbidden_patterns_from(other);
	}

	RMP& RMP::operator=(const RMP& other)
	{
		if (this != &other) {
			model = other.model;
			highs.passModel(model);
			copy_patterns_from(other);
			copy_forbidden_patterns_from(other);
		}
		return *this;
	}


	void RMP::copy_patterns_from(const RMP& other) {
		// 기존 패턴 삭제
		for (auto p : patterns)
			delete[] p;
		patterns.clear();

		// 새 패턴 복사
		for (auto src : other.patterns) {
			int* dst = new int[ProblemData::nL];
			for (int i = 0; i < ProblemData::nL; ++i)
				dst[i] = src[i];
			patterns.push_back(dst);
		}
	}

	void RMP::copy_forbidden_patterns_from(const RMP& other) {
		// 기존 forbidden pattern 정리 (혹시 잔재 있을 수 있음)
		for (auto pat : forbidden_patterns)
			delete[] pat;
		forbidden_patterns.clear();

		for (auto pat : other.forbidden_patterns) {
			int* new_pat = new int[ProblemData::nL];
			for (int i = 0; i < ProblemData::nL; ++i)
				new_pat[i] = pat[i];
			forbidden_patterns.push_back(new_pat);
		}
	}


	void RMP::initialize()
	{
		int n = ProblemData::nL;
		int L = ProblemData::LL;
		int* ReqL = ProblemData::ReqL;
		int* b = ProblemData::b;

		// initial patterns: trivial pattern for each order

		for (int i = 0; i < n; i++)
		{
			int* pat = new int[n]();
			pat[i] = floor(L / ReqL[i]);

			patterns.push_back(pat);
		}

		int num_vars = n;
		int num_cons = n;

		model.lp_.num_col_ = num_vars;
		model.lp_.num_row_ = num_cons;

		model.lp_.col_cost_.assign(num_vars, 1.0);
		model.lp_.col_lower_.assign(num_vars, 0.0);
		model.lp_.col_upper_.assign(num_vars, 1.0e30);

		model.lp_.row_lower_.assign(b, b + n);
		model.lp_.row_upper_.assign(num_cons, 1.0e30);

		vector<int> Astart(num_vars + 1);
		vector<int> Aindex;
		vector<double> Avalue;

		int nnz = 0;
		for (int p = 0; p < num_vars; ++p)
		{
			Astart[p] = nnz;
			for (int i = 0; i < n; ++i)
			{
				if (patterns[p][i] > 0)
				{
					Aindex.push_back(i);
					Avalue.push_back(patterns[p][i]);
					nnz++;
				}
			}
		}
		Astart[num_vars] = nnz;

		model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
		model.lp_.a_matrix_.start_ = Astart;
		model.lp_.a_matrix_.index_ = Aindex;
		model.lp_.a_matrix_.value_ = Avalue;

		highs.passModel(model);
	}

	double RMP::solve_RMP(vector<double>& duals)
	{
		for (auto ele : model.lp_.col_lower_)
		{
			cout << ele << "\t";
		}
		cout << endl;
		for (auto ele : model.lp_.col_upper_)
		{
			cout << ele << "\t";
		}
		cout << endl;

		highs.setOptionValue("output_flag", false);
		HighsStatus status = highs.run();
		HighsModelStatus model_status = highs.getModelStatus();

		if (model_status == HighsModelStatus::kInfeasible)
		{
			std::cout << "[RMP] Infeasible model detected." << std::endl;
			duals.assign(ProblemData::nL, 0.0);
			return std::numeric_limits<double>::infinity();
		}

		if (status != HighsStatus::kOk)
		{
			cerr << "HiGHS failed to solve RMP" << endl;
			return numeric_limits<double>::infinity();
		}

		HighsSolution sol = highs.getSolution();
		HighsInfo info = highs.getInfo();

		duals = sol.row_dual;
		double obj = info.objective_function_value;
		
		/*cout << "[RMP] LP objective: " << obj << endl;
		cout << "[RMP] duals: " << endl;
		for (int i = 0; i < duals.size(); ++i)
		{
			cout << duals[i] << " ";
		}*/

		return obj;
	}

	SP::SP()
	{
		SP_model.lp_.sense_ = ObjSense::kMaximize;
		highs.setOptionValue("output_flag", false);
	}

	bool SP::solve_SP(RMP& rmp, const vector<double>& duals)
	{
		int n = ProblemData::nL;
		int L = ProblemData::LL;
		int* ReqL = ProblemData::ReqL;
		int* b = ProblemData::b;

		SP_model.lp_.num_col_ = n;
		SP_model.lp_.num_row_ = 1;
		SP_model.lp_.col_cost_.assign(duals.begin(), duals.end());

		SP_model.lp_.col_lower_.assign(n, 0.0);
		SP_model.lp_.col_upper_.assign(n, 1.0e30);

		SP_model.lp_.row_lower_.push_back(-1.0e30);
		SP_model.lp_.row_upper_.push_back(L);

		vector<int> Astart(n + 1);
		vector<int> Aindex;
		vector<double> Avalue;

		Astart[0] = 0;
		for (int i = 0; i < n; ++i)
		{
			Aindex.push_back(0);
			Avalue.push_back(ReqL[i]);
			Astart[i + 1] = Aindex.size();
		}

		SP_model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
		SP_model.lp_.a_matrix_.start_ = Astart;
		SP_model.lp_.a_matrix_.index_ = Aindex;
		SP_model.lp_.a_matrix_.value_ = Avalue;

		SP_model.lp_.integrality_.assign(n, HighsVarType::kInteger);

		highs.setOptionValue("output_flag", false);
		highs.passModel(SP_model);
		HighsStatus status = highs.run();

		if (status != HighsStatus::kOk)
		{
			cerr << "Highs failed to solve SP" << endl;
			return true;
		}

		HighsSolution solution = highs.getSolution();
		HighsInfo info = highs.getInfo();

		vector<double> sol = solution.col_value;
		double obj = info.objective_function_value;

		if (1 - obj >= 0)
		{
			return true;
		}

		int* pattern = new int[n];

		for (int i = 0; i < n; ++i)
		{
			pattern[i] = static_cast<int>(sol[i] + 0.5);
		}

		if (is_forbidden_pattern(pattern, rmp.forbidden_patterns, n))
		{
			std::cout << "[sp] forbidden pattern generated. discarded." << std::endl;
			delete[] pattern;
			return true;
		}

		rmp.model.lp_.col_cost_.push_back(1.0);
		rmp.model.lp_.col_lower_.push_back(0.0);
		rmp.model.lp_.col_upper_.push_back(1.0e30);

		for (int i = 0; i < n; ++i)
		{
			if (pattern[i] > 0)
			{
				rmp.model.lp_.a_matrix_.index_.push_back(i);
				rmp.model.lp_.a_matrix_.value_.push_back(pattern[i]);
			}
		}

		rmp.patterns.push_back(pattern);
		rmp.model.lp_.a_matrix_.start_.push_back(rmp.model.lp_.a_matrix_.index_.size());
		rmp.model.lp_.num_col_ = rmp.model.lp_.col_cost_.size();
		rmp.highs.passModel(rmp.model);

		return false;
	}


	bool has_fractional_solution(RMP& rmp, int& branch_var, double& branch_value)
	{
		HighsSolution solution = rmp.highs.getSolution();
		vector<double> sol = solution.col_value;

		for (int i = 0; i < sol.size(); ++i)
		{
			double value = sol[i];
			if (abs(value - round(value)) > 1e-6)
			{
				branch_var = i;
				branch_value = value;
				return true;
			}
		}
		return false;
	}

	bool is_forbidden_pattern(int* new_pat, const vector<int*>& forbidden_patterns, int n)
	{
		for (auto pat : forbidden_patterns)
		{
			bool same = true;
			for (int i = 0; i < n; ++i)
			{
				if (new_pat[i] != pat[i])
				{
					same = false;
					break;
				}
			}
			if (same) return true;
		}
		return false;
	}


}