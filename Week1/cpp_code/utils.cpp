#include "utils.h"


// ------- STATISTICA -------

// calcolo della somma degli elementi di un vettore
double get_sum(vector<double>& v){
	double sum = 0;
	for(auto& el : v)
		sum += el;
	return sum;
}

// calcolo della media di un vettore
double get_mean(vector<double>& v){
	return get_sum(v)/v.size();
}

// calcolo della varianza di un vettore (media empirica)
double get_var(vector<double>& v){
	double sum = 0, m = get_mean(v);
	for(auto& el : v)
		sum += pow(el - m, 2);
	return sum/v.size();
}

// calcolo della varianza di un vettore (media teorica)
double get_var(vector<double>& v, double m){
	double sum = 0;
	for(auto& el : v)
		sum += pow(el - m, 2);
	return sum/v.size();
}


// calcolo della varianza di un vettore (media empirica)
double get_dev_std(vector<double>& v){
	return sqrt(get_var(v));
}

// calcolo della varianza di un vettore (media teorica)
double get_dev_std(vector<double>& v, double m){
	return sqrt(get_var(v, m));
}


// calcolo del chi quadro (su intervallo reale)
double get_chi_sq(vector<double>& v, double m, double min, double max, int n_int){
	
	double dx = (max - min)/double(n_int);
	vector<double> count;
	for(int i=0; i<n_int; i++)
		count.push_back(0);
	
	sort(v.begin(), v.end());
	for(auto& el : v){
		for(int i=0; i<n_int; i++){
         		if( el >= min + i*dx and el < min + (i+1)*dx ){
         			count[i] ++;
         			break;
         		}
         	}
	}
	
	return get_var(count, m) * n_int/m;
}

// salva un vettore
void save(vector<double>& v, string filename, string head){

	ofstream fout(filename);
	
	cout << "Inizio salvataggio su file: " << filename << endl;
	
	fout << head << endl;
	for(auto& el : v)
		fout << el << endl;
	
	fout.close();
	
	cout << "Fine salvataggio." << endl;
}



// ------- OPERAZIONI TRA VETTORI -------

//somma tra vettori
vector<double> operator+ (const vector<double>& v1, const vector<double>& v2){

	if(v1.size() != v2.size()){
		exit(0);
	}

	vector<double> v;
	for(int i=0; i<int(v1.size()); i++)
		v.push_back(v1[i] + v2[i]);

	return v;
	
}

//prodotto per costante
vector<double> operator* (double k, const vector<double>& v1){

	vector<double> v;
	for(auto& el : v1)
		v.push_back(k * el);

	return v;
	
}

//prodotto tra vettori
double operator* (const vector<double>& v1, const vector<double>& v2){

	if(v1.size() != v2.size()){
		exit(0);
	}

	double prod = 0;
	for(int i=0; i<int(v1.size()); i++)
		prod += v1[i] * v2[i];

	return prod;
	
}
