/**
 * @file algoritmoGeneroFijo.cpp
 * @brief Encuentra semigrupos numéricos internos y hojas para un género dado.
 * @author Mario Casas Pérez
 * @date 2025-06-16
 * @license MIT
 * @note Compilar con -std=c++11 o superior.
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <set>
#include <string>
#include <utility>
#include <chrono>
#include <regex>

using namespace std;

unordered_map<string, int> cacheGenero;
unordered_map<string, int> cacheFrobenius;

/**
 * @brief Calcula el máximo común divisor entre dos enteros.
 * @param[in] a Primer operando.
 * @param[in] b Segundo operando.
 * @return El máximo común divisor de a y b.
 */
int maxCD(int a, int b) {
   while (b) {
      a %= b;
      swap(a, b);
   }
   return a;
}

/**
 * @brief Genera una clave ordenada para memoización a partir de un vector de generadores.
 * @param[in] S Conjunto de generadores.
 * @return Cadena con los elementos de S ordenados y separados por comas.
 */
string clave(const vector<int>& S) {
   vector<int> copia = S;
   sort(copia.begin(), copia.end());
   string k = "";
   for (int x : copia) 
      k += to_string(x) + ",";
   return k;
}

/**
 * @brief Calcula el conductor de un semigrupo numérico.
 * @param[in] S Conjunto de generadores.
 * @return El menor n tal que todos los enteros >= n son representables.
 * @details
 *   Se expande el conjunto de alcanzables hasta encontrar un bloque de consecutivos
 *   de longitud igual al mínimo del conjunto de generadores.
 */
int calculaConductor(const vector<int>& S) {
   int minS = *min_element(S.begin(), S.end());
   int maxVal = accumulate(S.begin(), S.end(), 0);
   vector<bool> alcanzable;
   bool encontradoBloqueConsecutivo = false;
   int consecutivosNecesarios = minS;
   int conductor = maxVal;
   while (!encontradoBloqueConsecutivo) {
      alcanzable.assign(maxVal + 1, false);
      alcanzable[0] = true;

      for (int s : S)
         for (int i = s; i <= maxVal; ++i)
            alcanzable[i] = alcanzable[i] || alcanzable[i - s];

      int consecutivos = 0;
      for (int i = 0; i <= maxVal; ++i) {
         consecutivos = alcanzable[i] ? consecutivos + 1 : 0;
         if (consecutivos >= consecutivosNecesarios) {
            conductor = i - consecutivosNecesarios + 1;
            encontradoBloqueConsecutivo = true;
            break;
         }
      }
      if (!encontradoBloqueConsecutivo)
         maxVal *= 2;
   }
   return conductor;
}

/**
 * @brief Calcula el género de un semigrupo numérico.
 * @param[in] S Conjunto de generadores.
 * @return Número de enteros no representables por S.
 * @details
 *   Utiliza algoritmo similar al del conductor para determinar hasta dónde
 *   hay un bloque de consecutivos, y cuenta cuántos valores antes de ese
 *   conductor no son alcanzables.
 */
int calculaGenero(const vector<int>& S) {
   string k = clave(S);
   if (cacheGenero.count(k))
      return cacheGenero[k];

   int minS = *min_element(S.begin(), S.end());
   int maxVal = accumulate(S.begin(), S.end(), 0);
   vector<bool> alcanzable;

   bool encontradoBloqueConsecutivo = false;
   int consecutivosNecesarios = minS;
   int genero = 0;
   while (!encontradoBloqueConsecutivo) {
      alcanzable.assign(maxVal + 1, false);
      alcanzable[0] = true;

      for (int s : S)
         for (int i = s; i <= maxVal; ++i)
            alcanzable[i] = alcanzable[i] || alcanzable[i - s];

      int consecutivos = 0;
      for (int i = 0; i <= maxVal; ++i) {
         consecutivos = alcanzable[i] ? consecutivos + 1 : 0;
         if (consecutivos >= consecutivosNecesarios) {
            encontradoBloqueConsecutivo = true;
            break;
         }
      }
      if (!encontradoBloqueConsecutivo)
         maxVal *= 2;
   }

   int conductor = 0;
   {
      int consecutivos = 0;
      for (int i = 0; i <= maxVal; ++i) {
         consecutivos = alcanzable[i] ? consecutivos + 1 : 0;
         if (consecutivos >= consecutivosNecesarios) {
            conductor = i - consecutivosNecesarios + 1;
            break;
         }
      }
   }

   for (int i = 0; i < conductor; ++i) {
      if (!alcanzable[i])
         genero++;
   }
   return cacheGenero[k] = genero;
}

/**
 * @brief Calcula el número de Frobenius de un semigrupo generado por S.
 * @param[in] S Conjunto de generadores.
 * @return Máximo entero no representable (o -1 si todos son representables).
 */
int calculaFrobenius(const vector<int>& S) {
   string k = clave(S);
   if (cacheFrobenius.count(k))
      return cacheFrobenius[k];

   int maxVal = accumulate(S.begin(), S.end(), 0);
   vector<bool> alcanzable(maxVal + 1, false);
   alcanzable[0] = true;

   for (int s : S)
      for (int i = s; i <= maxVal; ++i)
         alcanzable[i] = alcanzable[i] || alcanzable[i - s];

   for (int i = maxVal; i >= 0; --i)
      if (!alcanzable[i])
         return cacheFrobenius[k] = i;

   return cacheFrobenius[k] = -1; // No existe número de Frobenius.
}

/**
 * @brief Comprueba si S es minimal en forma expandida de Hilbert.
 * @param[in] S Conjunto de generadores.
 * @return true si ningún generador es redundante y se alcanza un bloque mínimo.
 */
bool esMinimalHilbert(const vector<int>& S) {
   if (S.size() == 1) return true;

   int minS = *min_element(S.begin(), S.end());
   int conductor = calculaConductor(S);
   int limit = conductor + minS;

   set<int> hilbert;
   vector<int> cola = {0};
   hilbert.insert(0);

   while (!cola.empty()) {
      int t = cola.back(); 
      cola.pop_back();
      for (int a : S) {
         int nuevoTerm = t + a;
         if (nuevoTerm <= limit && !hilbert.count(nuevoTerm)) {
            hilbert.insert(nuevoTerm);
            cola.push_back(nuevoTerm);
         }
      }
   }

   int consecutivos = 0;
   int prev = -2; //el primero no se debe contar como consecutivo.
   for (int n : hilbert) {
      consecutivos = (n == prev + 1) ? consecutivos + 1 : 1;
      if (consecutivos >= minS) 
         break;
      prev = n;
   }
   if (consecutivos < minS) 
      return false;

   //Se comprueba que ningún elemento de S puede generarse sin él mismo.
   for (size_t i = 0; i < S.size(); ++i) {
      vector<int> copia = S;
      copia.erase(copia.begin() + i);
      set<int> formaExpandida;
      vector<int> cola = {0};
      formaExpandida.insert(0);
      int lim = S[i] + 1;
      while (!cola.empty()) {
         int t = cola.back(); 
         cola.pop_back();
         for (int a : copia) {
            int nuevo = t + a;
            if (nuevo <= lim && !formaExpandida.count(nuevo)) {
               formaExpandida.insert(nuevo);
               cola.push_back(nuevo);
            }
         }
      }
      if (formaExpandida.count(S[i]))
         return false; //El elemento puede generarse sin el mismo.
   }

   //Se comprueba que todos los elementos generadores están en la forma expandida.
   for (int s : S) {
      if (!hilbert.count(s))
         return false;
   }

   return true;
}

/**
 * @brief Comprueba si S define una hoja
 * @param[in] S Conjunto de generadores.
 * @param[in] genero Entero que indica el género.
 * @return true si S es hoja, false si no.
 */
bool esHoja(const vector<int>& S, int genero) {
   if (calculaGenero(S) != genero)
      return false;
   int f = calculaFrobenius(S);
   if (f == -1)
      return false;
   int max_elem = *max_element(S.begin(), S.end());
   return f > max_elem;
}

/**
 * @brief Comprueba si todos los pares en S son coprimos (mcd = 1).
 * @param[in] S Conjunto de generadores.
 * @return true si mcd de cada par es 1.
 */
bool mcdEsUno(const vector<int>& S) {
   if (S.empty()) return false;
   int mcd = S[0];
   for (size_t i = 1; i < S.size(); ++i) {
      mcd = maxCD(mcd, S[i]);
      if (mcd == 1) 
         return true;
   }
   return mcd == 1;
}

/**
 * @brief Genera todas las combinaciones de tamaño fijo a partir de un conjunto de números.
 * @param[out] subconjuntos Contenedor donde se almacenarán las combinaciones resultantes.
 * @param[in,out] combinacion Vector temporal que guarda la combinación en construcción.
 * @param[in] numeros Vector de entrada con los elementos de los que se extraen las combinaciones.
 * @param[in] comienzo Índice en `numeros` desde el cual se empieza a tomar elementos.
 * @param[in] tamano Tamaño objetivo de cada combinación.
 * @details
 * La función implementa un algoritmo de backtracking:
 *   1. Si la longitud de `combinacion` alcanza `tamano`, se añade una copia a `subconjuntos`.
 *   2. En caso contrario, recorre los índices desde `comienzo` hasta el final de `numeros`,
 *      añade cada elemento a `combinacion`, y llama recursivamente avanzando `comienzo`.
 *   3. Tras la llamada recursiva, se elimina el último elemento (backtrack) para explorar
 *      otras posibilidades.
 */
void generaCombinaciones(vector<vector<int>>& subconjuntos, vector<int>& combinacion, const vector<int>& numeros, int comienzo, int tamano) {
   if (combinacion.size() == tamano) {
      subconjuntos.push_back(combinacion);
      return;
   }
   for (int i = comienzo; i < numeros.size(); ++i) {
      combinacion.push_back(numeros[i]);
      generaCombinaciones(subconjuntos, combinacion, numeros, i + 1, tamano);
      combinacion.pop_back();
   }
}

/**
 * @brief Compara la cantidad de semigrupos numéricos internos y hojas, e imprime el resultado.
 * @param[in] internos Vector de semigrupos numéricos internos.
 * @param[in] hojas   Vector de semigrupos numéricos que son hojas.
 * @details
 * La función compara los tamaños de los dos contenedores:
 *   - Si hay más internos que hojas, informa que predominan los internos.
 *   - Si hay más hojas que internos, informa que predominan las hojas.
 *   - Si ambos tienen igual tamaño, informa que hay la misma cantidad.
 */
void comparaCantidades(const vector<vector<int>>& internos, const vector<vector<int>>& hojas) {
   cout << "\nComparacion:\n";
   if (internos.size() > hojas.size())
      cout << "Hay mas semigrupos numericos internos que hojas.\n";
   else if (internos.size() < hojas.size())
      cout << "Hay mas semigrupos numericos hojas que internos.\n";
   else
      cout << "Hay igual cantidad de semigrupos numericos internos y hojas.\n";
}

/**
 * @brief Encuentra e imprime todos los semigrupos internos y hojas para un género fijo.
 * @param genero Género fijo dado.
 */
void encontrarSemigruposYHojas(int genero) {
   vector<vector<int>> internos;
   vector<vector<int>> hojas;

   //se limita el espacio de búsqueda con límite.
   int limite = genero * 5;
   vector<int> numeros(limite - 1);
   iota(numeros.begin(), numeros.end(), 2);

   for (int tamano = 2; tamano <= genero; ++tamano) {
      vector<vector<int>> subconjuntos;
      vector<int> combinacion;
      generaCombinaciones(subconjuntos, combinacion, numeros, 0, tamano);

      for (const auto& subconjunto : subconjuntos) {
         if (!mcdEsUno(subconjunto))
            continue;
         if (calculaGenero(subconjunto) != genero)
            continue;
         if (!esMinimalHilbert(subconjunto))
            continue;

         if (esHoja(subconjunto, genero))
            hojas.push_back(subconjunto);
         else
            internos.push_back(subconjunto);
      }
   }
    
   vector<int> semigrupoExtra;
   for (int i = genero + 1; i <= 2 * genero + 1; ++i)
      semigrupoExtra.push_back(i);
   internos.push_back(semigrupoExtra);

   cout << "Semigrupos numericos internos:\n";
   for (const auto& s : internos) {
      cout << "<";
      for (size_t i = 0; i < s.size(); ++i)
         cout << s[i] << (i < s.size() - 1 ? "," : "");
      cout << ">\n";
   }

   cout << "\nSemigrupos numericos hoja:\n";
   for (const auto& h : hojas) {
      cout << "<";
      for (size_t i = 0; i < h.size(); ++i)
         cout << h[i] << (i < h.size() - 1 ? "," : "");
      cout << ">\n";
   }
    
   comparaCantidades(internos, hojas);
}

/**
 * @brief Función principal: lee el género dado por el usuario, realiza los cálculos e imprime los resultados y el tiempo.
 * @return 0 si éxito, 1 si error en entrada.
 */
int main() {
   int genero;
   cout << "Introduce el genero: ";
    
   string input;
   getline(cin, input);

   //Se comrpueba que la entrada sean solo números.
   regex patron("^[0-9]+$");
   if (!regex_match(input, patron)) {
      cout << "Entrada no valida. Debes introducir un unico numero entero (0 o mayor) sin caracteres." << endl;
      return 1;
   }

   //Cadena a número
   try {
      genero = stoi(input);
   } catch (...) {
      cout << "No se pudo convertir la entrada. Asegúrate de ingresar un número válido." << endl;
      return 1;
   }
    
   //Condición por defecto si el usuario introduce 0 como género.
   if (genero == 0) {
      cout << "Semigrupos numericos internos:\n";
      cout << "<1>\n";
      cout << "\nSemigrupos numericos hoja:\n";
      // No se muestran hojas
      return 0;
   }

   cout << "Calculando semigrupos numericos internos y hojas...\n";
    
   auto inicio = chrono::high_resolution_clock::now();
   encontrarSemigruposYHojas(genero);
   auto fin = chrono::high_resolution_clock::now();
    
   auto duracion = chrono::duration_cast<chrono::seconds>(fin - inicio).count();

   cout << "\nEl programa tardo " << duracion << " segundos.\n";

   return 0;
}

