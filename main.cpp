#include "BnP.h"
#include "ProblemData.h"
#include "Highs.h"
#include <iostream>
#include <stack>
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

int main()
{
	double best_obj = 1.0e30;
	vector<int*> best_patterns;
	vector<double> best_sol;
	
	bnp::RMP root;
	root.initialize();
	
	vector<double> duals;
	root.solve_RMP(duals);

	stack<bnp::RMP> node_stack;
	node_stack.push(move(root));

	int node_id = 0; 

	while (!node_stack.empty())
	{
		bnp::RMP current = move(node_stack.top());
		node_stack.pop();

		std::cout << "\n=============================================" << std::endl;
		std::cout << "[Node " << node_id++ << "] Processing new node." << std::endl;
		std::cout << "=============================================" << std::endl;

		//double current_obj;
		double current_obj = current.solve_RMP(duals);

		// Column Generation
		bool col_gen_term = false;
		bnp::SP sp;
		while (true)
		{
			col_gen_term = sp.solve_SP(current, duals);
			if (col_gen_term)
			{
				std::cout << "[Column Generation] No improving pattern found. Terminated." << std::endl;
				break;
			}

			current_obj = current.solve_RMP(duals);
		}
		/*if (col_gen_term)
			current_obj = current.solve_RMP(duals);*/

		std::cout << "\n[LP BOUND] Objective value (relaxation): " << current_obj << std::endl;

		if (current_obj >= best_obj)
		{
			std::cout << "[Prune] Node pruned by bound (current_obj >= best_obj)." << std::endl;
			continue;
		}

		std::cout << "\n[Patterns in this node] Total: " << current.patterns.size() << endl;
		for (int p = 0; p < current.patterns.size(); ++p)
		{
			std::cout << "Pattern " << p << ": [";
			for (int i = 0; i < ProblemData::nL; ++i)
			{
				std::cout << current.patterns[p][i] << " ";
			}
			std::cout << "]" << std::endl;
		}

		auto current_sol = current.highs.getSolution().col_value;
		std::cout << "\n[LP Solution] Pattern usages:" << std::endl;
		for (int i = 0; i < current_sol.size(); ++i)
		{
			std::cout << "x_" << i << " = " << current_sol[i] << std::endl;
		}

		// check if the solution is integral
		int branch_var;
		double branch_value;
		bool fractional = bnp::has_fractional_solution(current, branch_var, branch_value);
		
		std::cout << "\n[Solution Status] Is integer feasible? " << (fractional ? "No" : "Yes") << std::endl;

		if (fractional)
		{
			std::cout << "-> Fractional variable: x_" << branch_var << " = " << branch_value << std::endl;
		}

		if (fractional && (current_obj < best_obj))
		{
			// branching
			bnp::RMP left(current);
			left.model.lp_.col_lower_[branch_var] = 0.0;
			left.model.lp_.col_upper_[branch_var] = floor(branch_value);
			/*left.model.lp_.col_lower_[branch_var] = std::max(left.model.lp_.col_lower_[branch_var], 0.0);
			left.model.lp_.col_upper_[branch_var] = std::min(left.model.lp_.col_upper_[branch_var], floor(branch_value));*/

			left.highs.passModel(left.model);

			int* forbidden = new int[ProblemData::nL];
			for (int i = 0; i < ProblemData::nL; ++i)
				forbidden[i] = current.patterns[branch_var][i];
			left.forbidden_patterns.push_back(forbidden);

			node_stack.push(move(left));

			bnp::RMP right(current);
			right.model.lp_.col_lower_[branch_var] = ceil(branch_value);
			right.model.lp_.col_upper_[branch_var] = 1.0e30;
			/*right.model.lp_.col_lower_[branch_var] = std::max(right.model.lp_.col_lower_[branch_var], ceil(branch_value));
			right.model.lp_.col_upper_[branch_var] = std::min(right.model.lp_.col_upper_[branch_var], 1.0e30);*/
			right.highs.passModel(right.model);
			right.forbidden_patterns = {};
			
			node_stack.push(move(right));

		}

		if (!fractional && (current_obj < best_obj))
		{
			best_obj = current_obj;
			//best_patterns = current.patterns;

			best_patterns.clear();
			for (auto pat : current.patterns)
			{
				int* new_pat = new int[ProblemData::nL];
				for (int i = 0; i < ProblemData::nL; ++i)
					new_pat[i] = pat[i];
				best_patterns.push_back(new_pat);
			}

			best_sol = current.highs.getSolution().col_value;

			cout << "new incumbent found: " << best_obj << endl;
		}
	}

	// print result
	cout << "\n==============================\n";
	cout << "minimum number of rolls: " << best_obj << endl;
	
	for (int p = 0; p < best_patterns.size(); ++p)
	{
		if (best_sol[p] < 1e-6) continue;
		cout << "Pattern " << p << " (" << best_sol[p] << " times used) : [";
		for (int i = 0; i < ProblemData::nL; ++i)
		{
			cout << best_patterns[p][i] << " ";
		}
		cout << "]";
		cout << endl;
	}

	for (auto pat : best_patterns)
		delete[] pat;

	return 0;
}