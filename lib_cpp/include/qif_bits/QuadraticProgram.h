namespace qp {

using std::string;

// Use a map <col, row> => value for sparse arrays. Put column first so that
// the lexicographic order used by map is top-to-bottom/left-to-right
template<typename eT> using Sparse = std::map<std::pair<uint,uint>,eT>;

enum class Status { OPTIMAL, INFEASIBLE, ERROR };
enum class Method { ADDM };

std::ostream& operator<<(std::ostream& os, const Status& status);
std::ostream& operator<<(std::ostream& os, const Method& method);


class Defaults {
	public:
		static OSQPInt osqp_polish;
		static OSQPInt osqp_verbose;
		static OSQPFloat osqp_alpha;
		static OSQPFloat osqp_eps_abs;
		static OSQPFloat osqp_eps_rel;
		static Method	method;
};

// Solve the quadratic program
// min 1/2 x^T P x + c^T x
// subject to l <= A x <= u
// x >= 0

template<typename eT>
class QuadraticProgram {
	typedef uint Var;
	typedef uint Con;

	public:
		bool non_negative = false;
		Method method = Defaults::method;
		Status status;

		OSQPInt osqp_polish = Defaults::osqp_polish;
		OSQPInt osqp_verbose = Defaults::osqp_verbose;
		OSQPFloat osqp_alpha = Defaults::osqp_alpha;
		OSQPFloat osqp_eps_abs = Defaults::osqp_eps_abs;
		OSQPFloat osqp_eps_rel = Defaults::osqp_eps_rel;

		QuadraticProgram() {}

		bool solve();

		inline eT objective()				{ return obj; }
		inline eT solution(Var x)			{ return sol(x); }
		inline Col<eT> solution()			{ return sol; };

		void clear();
		void from_matrix(const arma::SpMat<eT>&P, const Col<eT>& c, const arma::SpMat<eT>& A, const Col<eT>& l, const Col<eT>& u);

		// API for creating variables and constraints
		Var make_var();
		std::vector<Var> make_vars(uint n);
		std::vector<std::vector<Var>> make_vars(uint n, uint m);

		Con make_con(eT lb, eT ub);

		void set_obj_coeff(Var var, eT coeff, bool add = false);				// linear part
		void set_obj_coeff(Var var1, Var var2, eT coeff, bool add = false);		// quadratic part
		void set_con_coeff(Con cons, Var var, eT coeff, bool add = false);

	protected:
		Col<eT> sol;			// solution
		eT obj;					// objective

		uint n_var = 0,							// number of variables
			 n_con = 0;							// number of constraints

		std::vector<OSQPFloat> obj_coeff_lin;		// coefficients for the objective function, linear part
		Sparse<eT> obj_coeff_quad;				// coefficients for the objective function, quadratic part
		Sparse<eT> con_coeff;					// coefficients for the constraints
		std::vector<OSQPFloat> con_lb, con_ub;	// constraints lower/upper

		bool osqp();
};

template<typename eT>
inline
typename QuadraticProgram<eT>::Var QuadraticProgram<eT>::make_var() {
	obj_coeff_lin.push_back(0);
	return n_var++;
}

template<typename eT>
inline
std::vector<typename QuadraticProgram<eT>::Var> QuadraticProgram<eT>::make_vars(uint n) {
	std::vector<Var> res(n);
	for(uint i = 0; i < n; i++)
		res[i] = make_var();
	return res;
}

template<typename eT>
inline
std::vector<std::vector<typename QuadraticProgram<eT>::Var>> QuadraticProgram<eT>::make_vars(uint n, uint m) {
	std::vector<std::vector<Var>> res(n);
	for(uint i = 0; i < n; i++) {
		res[i].resize(m);
		for(uint j = 0; j < m; j++)
			res[i][j] = make_var();
	}
	return res;
}

template<typename eT>
inline
typename QuadraticProgram<eT>::Con QuadraticProgram<eT>::make_con(eT lb, eT ub) {
	if(ub == infinity<eT>() && lb == -ub)
		throw std::runtime_error("trying to add unconstrained constraint");

	con_lb.push_back(lb);
	con_ub.push_back(ub);
	return n_con++;
}

template<typename eT>
inline
void QuadraticProgram<eT>::set_obj_coeff(Var var, eT coeff, bool add) {
	if(add)
		obj_coeff_lin[var] += coeff;
	else
		obj_coeff_lin[var] = coeff;
}

template<typename eT>
inline
void QuadraticProgram<eT>::set_obj_coeff(Var var1, Var var2, eT coeff, bool add) {
	if(equal<eT>(coeff, eT(0)))
		return;

	// maintain upper-triangular
	if(var1 > var2)
		std::swap(var1, var2);

	auto key = std::pair(var2, var1);		// <col, row>
	if(add && obj_coeff_quad.count(key))
		obj_coeff_quad[key] += coeff;
	else
		obj_coeff_quad[key] = coeff;
}

template<typename eT>
inline
void QuadraticProgram<eT>::set_con_coeff(Con con, Var var, eT coeff, bool add) {
	if(equal<eT>(coeff, eT(0)))
		return;

	auto key = std::pair(var, con);	// <col, row>
	if(add && con_coeff.count(key))
		con_coeff[key] += coeff;
	else
		con_coeff[key] = coeff;
}

template<typename eT>
void QuadraticProgram<eT>::clear() {
	obj_coeff_lin.clear();
	obj_coeff_quad.clear();
	con_coeff.clear();
	con_lb.clear();
	con_ub.clear();
	n_var = n_con = 0;
}

// input program in matrix form
// P: objective coeff, quadratic part
// c: objective constants, linear part
// A: constraint coefficients
// l: lower bounds for constraints
// u: upper bounds for constraints
//
template<typename eT>
void QuadraticProgram<eT>::from_matrix(const arma::SpMat<eT>&P, const Col<eT>& c, const arma::SpMat<eT>& A, const Col<eT>& l, const Col<eT>& u) {

	clear();
	n_var = A.n_cols;
	n_con = A.n_rows;

	if(P.n_rows != n_var || P.n_cols != n_var || c.n_elem != n_var || l.n_elem != n_con || u.n_elem != n_con)
		throw std::runtime_error("invalid size");
	
	auto end = P.end();
	for(auto c = P.begin(); c != end; ++c) {		// c++ throws weird warning, ++c doesn't!
		set_obj_coeff(c.row(), c.col(), *c);
	}

	using cfv = std::vector<OSQPFloat>;
	obj_coeff_lin = arma::conv_to<cfv>::from(c);
	
	end = A.end();
	for(auto c = A.begin(); c != end; ++c) {		// c++ throws weird warning, ++c doesn't!
		set_con_coeff(c.row(), c.col(), *c);
	}

	con_lb = arma::conv_to<cfv>::from(l);
	con_ub = arma::conv_to<cfv>::from(u);
}

template<typename eT>
bool QuadraticProgram<eT>::solve() {
	return osqp();
}

template<typename eT>
OSQPCscMatrix* to_csc(uint n_rows, uint n_cols, const Sparse<eT>& entries) {
	// Compressed Sparse Column (CSC) format.
	// https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
	// val:      array of non-zero values, in top-to-bottom, left-to-right order
	// row_ind:  for every value, the row index of that value
	// col_ptr:  for every column, the val index of the first element of that column.
	//           Plus an extra element at the end, with the size of val

	uint n_nonzero = entries.size();

	OSQPFloat* val = (OSQPFloat*) malloc(sizeof(OSQPFloat) * n_nonzero);
	OSQPInt* row_ind = (OSQPInt*) malloc(sizeof(OSQPInt) * n_nonzero);
	OSQPInt* col_ptr = (OSQPInt*) malloc(sizeof(OSQPInt) * n_cols + 1);

	uint cur_col = 0;
	uint cnt = 0;		// next to update
	for(auto& [key, value] : entries) {
		auto& [col, row] = key;

		// first time we see a column, update col_ptr
		for(; cur_col <= col; cur_col++)	// possibly update previous empty columns
			col_ptr[cur_col] = cnt;

		val[cnt] = value;
		row_ind[cnt] = row;
		cnt++;
	}
	for(; cur_col <= n_cols; cur_col++)
		col_ptr[cur_col] = n_nonzero;

	return wrapper::OSQPCscMatrix_new(n_rows, n_cols, n_nonzero, val, row_ind, col_ptr);
}


template<typename eT>
bool QuadraticProgram<eT>::osqp() {

	// TODO
	if(non_negative)
		throw std::runtime_error("not_implemented");

	// populate data
	//OSQPData* data = (OSQPData *)malloc(sizeof(OSQPData));
	//data->n = n_var;								// number of variables
	//data->m = n_con;								// number of constraints
	//data->q = obj_coeff_lin.data();					// cost function, linear part
	//data->P = to_csc(n_var, n_var, obj_coeff_quad);	// cost function, quadratic part
	//data->A = to_csc(n_con, n_var, con_coeff);		// constraints
	//data->l = con_lb.data();						// lower bounds
	//data->u = con_ub.data();						// upper bounds

	// settings
	OSQPSettings* settings = OSQPSettings_new();
	wrapper::osqp_set_default_settings(settings);
	settings->polishing = osqp_polish;
	settings->verbose = osqp_verbose;
	if(osqp_alpha   >= 0.0) settings->alpha = osqp_alpha;
	if(osqp_eps_abs >= 0.0) settings->eps_abs = osqp_eps_abs;
	if(osqp_eps_rel >= 0.0) settings->eps_rel = osqp_eps_rel;

	// solve
	OSQPSolver* solver;
  OSQPCscMatrix* P = to_csc(n_var, n_var, obj_coeff_quad);
  OSQPCscMatrix* A = to_csc(n_con, n_var, con_coeff);

	wrapper::osqp_setup(&solver, 
      P, // cost function, quadratic part
      obj_coeff_lin.data(), // cost function, linear part
      A, // constraints
      con_lb.data(), // lower bounds
      con_ub.data(), // upper bounds
      n_con, // number of constraints
      n_var, // number of variables
      settings);
	wrapper::osqp_solve(solver);
	// std::cout << "DONE\n";

	OSQPInt st = solver->info->status_val;
	status =
		st == OSQP_SOLVED											? Status::OPTIMAL :
		st == OSQP_PRIMAL_INFEASIBLE || st == OSQP_DUAL_INFEASIBLE	? Status::INFEASIBLE :
		Status::ERROR;

	// get optimal solution
	if(status == Status::OPTIMAL) {
		sol.set_size(n_var);
		for(uint j = 0; j < n_var; j++)
			sol.at(j) = solver->solution->x[j];

		obj = solver->info->obj_val;
	}

	// cleanup
  OSQPCscMatrix_free(P);
  OSQPCscMatrix_free(A);
	OSQPSettings_free(settings);
	wrapper::osqp_cleanup(solver);

	return status == Status::OPTIMAL;
}


} // namespace qp

