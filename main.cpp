#include "BnP.h"
#include "ProblemData.h"
#include "Highs.h"
#include <iostream>
#include <stack>
#include <vector>
#include <limits>

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

	while (!node_stack.empty())
	{
		bnp::RMP& current = move(node_stack.top());
		node_stack.pop();

		double current_obj;

		// Column Generation
		bool col_gen_term = false;
		bnp::SP sp;
		while (true)
		{
			col_gen_term = sp.solve_SP(current, duals);
			if (col_gen_term)
				break;

			current_obj = current.solve_RMP(duals);
		}

		// check if the solution is integral
		int branch_var;
		double branch_value;
		bool fractional = bnp::has_fractional_solution(current, branch_var, branch_value);
		
		if (fractional && (current_obj < best_obj))
		{
			// branching
			bnp::RMP left;
			left.model = current.model;
			left.patterns = current.patterns;
			left.model.lp_.col_lower_[branch_var] = 0.0;
			left.model.lp_.col_upper_[branch_var] = floor(branch_value);
			left.highs.passModel(left.model);
			node_stack.push(move(left));

			bnp::RMP right;
			right.model = current.model;
			right.patterns = current.patterns;
			right.model.lp_.col_lower_[branch_var] = ceil(branch_value);
			right.model.lp_.col_upper_[branch_var] = 1.0e30;
			right.highs.passModel(right.model);
			node_stack.push(move(right));
		}

		if (!fractional && (current_obj < best_obj))
		{
			best_obj = current_obj;
			best_patterns = current.patterns;

			best_sol = current.highs.getSolution().col_value;

			cout << "new incumbent found: " << best_obj << endl;
		}
	}

	// print result
	cout << "\n==============================\n";
	cout << "최소 롤 수: " << best_obj << endl;
	
	for (int p = 0; best_patterns.size(); ++p)
	{
		if (best_sol[p] < 1e-6) continue;
		cout << "패턴 " << p << " (" << best_sol[p] << "번 사용) : [";
		for (int i = 0; i < ProblemData::nL; ++i)
		{
			cout << best_patterns[p][i] << " ";
		}
		cout << endl;
	}
}