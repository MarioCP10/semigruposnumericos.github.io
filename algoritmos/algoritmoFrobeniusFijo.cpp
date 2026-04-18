/**
 * @file algoritmoFrobeniusFijo.cpp
 * @brief Genera e imprime todos los semigrupos numéricos con un valor fijo de Frobenius.
 * @details
 * El programa construye el semigrupo inicial S0 = <F+1, F+2, …, 2F+1>, lo minimiza,
 * y luego realiza una búsqueda en amplitud añadiendo de forma recursiva nuevos generadores
 * que mantengan el mismo Frobenius F. Para cada semigrupo encontrado, muestra su 
 * representación en forma de lista y su conjunto de Apéry respecto a F+1.
 * Finalmente se clasifica y lista cuáles son semigrupos numéricos internos y
 * cuáles son semigrupos numéricos hojas.
 * @author Mario Casas Pérez
 * @date 2025-06-16
 * @license MIT
 * @note Compilar con -std=c++11 o superior.
 */

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>  //Se incluye para formateo (setw, left).
#include <chrono>
#include <regex>

using namespace std;

/**
 * @brief Convierte un vector de generadores en su notación de semigrupo.
 * @param[in] generadores Conjunto de generadores del semigrupo.
 * @return Cadena con el semigrupo en notación angular.
 */
string semigrupoAString(const vector<int>& generadores) {
   ostringstream oss;
   oss << "<";
   for (size_t i = 0; i < generadores.size(); i++) {
      oss << generadores[i];
      if (i + 1 < generadores.size()) 
         oss << ",";
   }
   oss << ">";
   return oss.str();
}

/**
 * @brief Igual que semigrupoAString; alias usado para claridad en contextos de Apéry.
 * @param[in] generadores Conjunto de generadores.
 * @return Misma salida que semigrupoAString.
 */
string corchetesAngulares(const vector<int>& generadores) {
   ostringstream oss;
   oss << "<";
   for (size_t i = 0; i < generadores.size(); i++) {
      oss << generadores[i];
      if (i + 1 < generadores.size()) 
         oss << ",";
   }
   oss << ">";
   return oss.str();
}

/**
 * @brief Comprueba si un valor es representable como combinación no negativa de generadores.
 * @param[in] valor Entero a verificar.
 * @param[in] generadores Conjunto de generadores del semigrupo.
 * @return true si existe combinación de enteros >= 0 que sume `valor`; false en caso contrario.
 * @details
 * Utiliza programación dinámica: dp[i] = true si i es alcanzable. Inicializa dp[0]=true,
 * y para cada i alcanzable, marca dp[i + g] para cada g en generadores.
 */
bool esRepresentable(int valor, const vector<int>& generadores) {
   if (valor < 0) 
      return false;
   vector<bool> dp(valor + 1, false);
   dp[0] = true;
   for (int i = 0; i <= valor; i++) {
      if (!dp[i]) continue;
      for (int g : generadores) {
         if (i + g <= valor) 
            dp[i + g] = true;
      }
   }
   return dp[valor];
}

/**
 * @brief Reduce un conjunto de generadores eliminando los redundantes.
 * @param[in] generadores Vector inicial de generadores.
 * @return Nuevo vector donde se han eliminado aquellos g que sean representables
 *         por el resto de generadores.
 * @details
 * Ordena primero los generadores para procesarlos de menor a mayor. Para cada g,
 * lo elimina temporalmente y comprueba si g sigue siendo representable; si no lo es,
 * se conserva en el conjunto mínimo.
 */
vector<int> minimizarGeneradores(vector<int> generadores) {
   sort(generadores.begin(), generadores.end());
   vector<int> minimal;
   for (size_t i = 0; i < generadores.size(); i++) {
      vector<int> otros = generadores;
      otros.erase(otros.begin() + i);
      if (!esRepresentable(generadores[i], otros))
         minimal.push_back(generadores[i]);
   }
   return minimal;
}

/**
 * @brief Verifica la validez del Frobenius de un semigrupo.
 * @param[in] generadores Conjunto de generadores minimizados.
 * @param[in] F Valor de Frobenius a comprobar.
 * @return true si F no es representable pero F+1 sí lo es; false en otro caso.
 */
bool frobeniusValido(const vector<int>& generadores, int F) {
   return !esRepresentable(F, generadores)
      &&  esRepresentable(F+1, generadores);
}

/**
 * @brief Construye el semigrupo inicial S0 = <F+1, F+2, …, 2F+1>.
 * @param[in] F Valor de Frobenius.
 * @return Vector de generadores de S0.
 */
vector<int> semigrupoInicial(int F) {
   vector<int> generadores;
   for (int i = F+1; i <= 2*F+1; i++)
      generadores.push_back(i);
   return generadores;
}

/**
 * @brief Representa un candidato al agregar un nuevo generador.
 */
struct Candidato {
   int x;   //Elemento candidato a añadir.
   vector<int> semigrupo;  //Semigrupo resultante tras la minimización.
};

/**
 * @brief Obtiene los candidatos que pueden añadirse manteniendo Frobenius = F.
 * @param[in] S Conjunto actual de generadores minimizados.
 * @param[in] F Valor de Frobenius objetivo.
 * @return Vector de estructuras Candidato con cada posible x y su semigrupo.
 */
vector<Candidato> obtenerCandidatos(const vector<int>& S, int F) {
   vector<Candidato> candidatos;
   int m = S.empty() ? F+1 : S[0];
   for (int x = 2; x < m; x++) {
      if (x == F) 
         continue;
      vector<int> T = S;
      if (find(T.begin(), T.end(), x) == T.end()) {
         T.push_back(x);
         T = minimizarGeneradores(T);
         if (frobeniusValido(T, F))
            candidatos.push_back({x, T});
      }
   }
   return candidatos;
}

/**
 * @brief Calcula el conjunto de Apéry de un semigrupo respecto a un período.
 * @param[in] generadores Conjunto de generadores minimizados.
 * @param[in] periodo     Entero de referencia para el resto (por lo general F+1).
 * @return Vector ordenado de size=periodo con el menor valor representable para cada resto.
 * @details
 * Emplea programación dinámica similar a esRepresentable hasta 2·periodo² para asegurar
 * que cada clase de resto mod periodo aparezca.
 */
vector<int> conjuntoApery(const vector<int>& generadores, int periodo) {
   int limite = 2 * periodo * periodo;
   vector<bool> dp(limite + 1, false);
   dp[0] = true;
   for (int i = 0; i <= limite; i++) {
      if (!dp[i]) continue;
      for (int g : generadores) {
         int nx = i + g;
         if (nx <= limite) 
            dp[nx] = true;
      }
   }
   vector<int> ap;
   for (int r = 0; r < periodo; r++) {
      int candidato = -1;
      for (int v = 0; v <= limite; v++) {
         if (dp[v] && v % periodo == r) {
            candidato = v;
            break;
         }
      }
      if (candidato >= 0) ap.push_back(candidato);
   }
   sort(ap.begin(), ap.end());
   return ap;
}

/**
 * @brief Convierte el conjunto de Apéry a cadena.
 * @param[in] generadores Conjunto de generadores.
 * @param[in] periodo     Período usado en el cálculo (F+1).
 * @return Cadena con formato Ap(<generadores>, periodo) = {a0,…,a_{periodo-1}}.
 */
string aperyAString(const vector<int>& generadores, int periodo) {
   auto ap = conjuntoApery(generadores, periodo);
   ostringstream oss;
   oss << "Ap(" << corchetesAngulares(generadores) << ", " << periodo << ") = {";
   for (size_t i = 0; i < ap.size(); i++) {
      oss << ap[i];
      if (i + 1 < ap.size()) oss << ",";
   }
   oss << "}";
   return oss.str();
}

/**
 * @brief Comprueba si un semigrupo ya ha sido procesado para evitar duplicados.
 * @param[in] S     Conjunto de generadores minimizados.
 * @param[in] vistos Conjunto de representaciones en cadena ya mostradas.
 * @return true si semigrupoAString(S) está en vistos.
 */
bool yaVisto(const vector<int>& S, const set<string>& vistos) {
   return vistos.count( semigrupoAString(S) ) > 0;
}

/**
 * @brief Genera e imprime todos los semigrupos con Frobenius = F en BFS.
 * @param[in] F Valor de Frobenius deseado.
 * @return Vector de todos los semigrupos minimizados encontrados.
 * @details
 * - Inicia con S0 minimizado.
 * - Repite: para cada semigrupo en el nivel actual, genera candidatos, filtra duplicados,
 *   los imprime alineados con su Apéry y los añade al siguiente nivel.
 */
vector<vector<int>> generaSemigruposConF(int F) {
   set<string> vistos;
   vector<vector<int>> resultado;

   //Semigrupo inicial
   vector<int> S0 = semigrupoInicial(F);
   S0 = minimizarGeneradores(S0);
   resultado.push_back(S0);
   vistos.insert(semigrupoAString(S0) );

   string s0_str = semigrupoAString(S0);
   string s0_apery = aperyAString(S0, F+1);
   cout << left << setw(30) << s0_str << " | " << s0_apery << "\n";

   vector<vector<int>> nivelActual{ S0 };

   //Se hace una búsqueda en amplitud sobre las ramas
   while (!nivelActual.empty()) {
      vector<vector<int>> siguienteNivel;
      for (auto &S : nivelActual) {
         auto candidatos = obtenerCandidatos(S, F);
         for (auto &c : candidatos) {
            if (!yaVisto(c.semigrupo, vistos)) {
               vistos.insert( semigrupoAString(c.semigrupo) );
               resultado.push_back(c.semigrupo);
               siguienteNivel.push_back(c.semigrupo);
               //Se muestra el semigrupo y el conjunto de Apéry en columnas alineadas
               string semigrupoStr = semigrupoAString(c.semigrupo);
               string aperyStr = aperyAString(c.semigrupo, F+1);
               cout << left << setw(30) << semigrupoStr << " | " << aperyStr << "\n";
            }
         }
      }
      nivelActual.swap(siguienteNivel);
   }

   return resultado;
}

/**
 * @brief Punto de entrada: valida F, lanza la generación y clasifica resultados.
 * @return Código de salida (0 éxito, 1 error de entrada).
 */
int main(){
   string input;
   cout << "Numero de Frobenius (F): ";
   getline(cin, input);

   //Se comrpueba que la entrada sean solo números (un único entero positivo)
   regex pattern("^[0-9]+$");
   if (!regex_match(input, pattern)) {
      cout << "Entrada no valida. Debes introducir un unico numero entero positivo (sin decimales, letras o espacios extras)." << endl;
      return 1;
   }

   int F = stoi(input);

   if (F <= 0) {
      cout << "F debe ser un entero positivo." << endl;
      return 1;
   }
    
   //Se generan todos los semigrupos numéricos con Frobenius = F
   auto inicio = chrono::high_resolution_clock::now();
    
   auto todos = generaSemigruposConF(F);
    
   //Se pasa a clasificar si los semigrupos numéricos son internos u hojas
   vector<vector<int>> internos, hojas;
   for (auto &S : todos) {
      bool todoMenor = true;
      for (int g : S) {
         if (g >= F) {
            todoMenor = false;
            break;
         }
      }
      if (todoMenor)
         hojas.push_back(S);
      else
         internos.push_back(S);
   }
    
   auto fin = chrono::high_resolution_clock::now();
   auto duracion = chrono::duration_cast<chrono::seconds>(fin - inicio).count();

   cout << "\nSemigrupos numericos internos\n";
   for (auto &S : internos) {
      cout << semigrupoAString(S) << "\n";
   }

   cout << "\nSemigrupos numericos hoja\n";
   for (auto &S : hojas) {
      cout << semigrupoAString(S) << "\n";
   }

   cout << "\nTotal internos: " << internos.size()
      << "   Total hojas: " << hojas.size() << "\n";
      
   cout << "\nEl programa ha tardado " << duracion << " segundos.\n";

   return 0;
}

