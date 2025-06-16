/**
 * @file algoritmoGeneroMultiplicidadFija.cpp
 * @brief Encuentra semigrupos numéricos internos y hojas para un género y multiplicidad dados.
 * @details
 *   El programa genera todas las combinaciones de generadores donde la multiplicidad
 *   es especificada, comprobando que generen el género deseado. Clasifica los semigrupos en
 *   internos y hojas, imprimiendo ambos conjuntos y comparando sus tamaños.
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
   if (S.size() == 1) 
      return true;

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
 * @brief Determina si S genera un semigrupo hoja (F > max(S)).
 * @param[in] S Conjunto de generadores.
 * @param[in] genero Entero que indica el género.
 * @return true si S genera el género dado y su Frobenius supera al máximo generador.
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
 * @brief Genera recursivamente todas las combinaciones de tamaño tamano.
 * @param[out] subconjuntos Acumula las combinaciones generadas.
 * @param[in,out] combinacion Vector auxiliar para la combinación actual.
 * @param[in] numeros Conjunto de números candidatos.
 * @param[in] comienzo Índice desde el que añadir en esta llamada.
 * @param[in] tamano Longitud deseada de cada combinación.
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
 * @brief Comprueba si los elementos de S son globalmente coprimos.
 * @param[in] S Conjunto de generadores.
 * @return true si el mcd acumulado termina en 1.
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
 * @brief Compara y muestra cuál de dos conjuntos tiene más semigrupos.
 * @param[in] internos Lista de semigrupos numéricos internos.
 * @param[in] hojas    Lista de semigrupos numéricos hoja.
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
 * @brief Busca y clasifica semigrupos numéricos de género y multiplicidad fijos.
 * @param[in] genero        Género objetivo.
 * @param[in] multiplicidad Mínimo elemento (multiplicidad) de cada semigrupo sin contar el 0.
 * @details
 *   - Genera todas las combinaciones de tamaño 2..género de los números
 *     desde multiplicidad hasta 5·género.
 *   - Filtra por mcdEsUno, cálculo exacto de género y minimalidad de Hilbert.
 *   - Si cumple esHoja, entonces hoja; si no, interno.
 *   - Si multiplicidad == género+1, añade el semigrupo <F+1,…,2F+1>.
 *   - Imprime ambos conjuntos y compara sus tamaños.
 */
void encontrarSemigruposYHojas(int genero, int multiplicidad) {
   vector<vector<int>> internos;
   vector<vector<int>> hojas;

   //se limita el espacio de búsqueda con límite
   int limite = genero * 5;
   vector<int> numeros;
   for (int i = multiplicidad; i <= limite; ++i)
      numeros.push_back(i);

   for (int tamano = 2; tamano <= genero; ++tamano) {
      vector<vector<int>> subconjuntos;
      vector<int> combinacion;
      generaCombinaciones(subconjuntos, combinacion, numeros, 0, tamano);

      for (const auto& subconjunto : subconjuntos) {
         //se debe incluir la multiplicidad como mínimo
         if (*min_element(subconjunto.begin(), subconjunto.end()) != multiplicidad)
            continue;
         if (!mcdEsUno(subconjunto))
            continue;
         if (calculaGenero(subconjunto) != genero)
            continue;
         if (!esMinimalHilbert(subconjunto))
            continue;

         //Se clasifica en semigrupos numéricos internos u hojas
         if (esHoja(subconjunto, genero))
            hojas.push_back(subconjunto);
         else
            internos.push_back(subconjunto);
      }
   }
    
   if(multiplicidad == genero + 1){
      vector<int> semigrupoExtra;
      for (int i = genero + 1; i <= 2 * genero + 1; ++i)
         semigrupoExtra.push_back(i);
      internos.push_back(semigrupoExtra);
   }

   cout << "\nSemigrupos numericos internos (m=" << multiplicidad << ", g=" << genero << "):\n";
   for (const auto& s : internos) {
      cout << "<";
      for (size_t i = 0; i < s.size(); ++i)
         cout << s[i] << (i + 1 < s.size() ? "," : "");
      cout << ">\n";
   }

   cout << "\nSemigrupos numericos hoja (m=" << multiplicidad << ", g=" << genero << "):\n";
   for (const auto& h : hojas) {
      cout << "<";
      for (size_t i = 0; i < h.size(); ++i)
         cout << h[i] << (i + 1 < h.size() ? "," : "");
      cout << ">\n";
   }

   comparaCantidades(internos, hojas);
}

/**
 * @brief Punto de entrada: valida parámetros y se ejecuta el programa.
 * @return 0 si éxito, 1 si error en entrada.
 */
int main() {
   string inputGenero, inputMultiplicidad;
    
   cout << "Introduce el genero: ";
   getline(cin, inputGenero);
    
   cout << "Introduce la multiplicidad: ";
   getline(cin, inputMultiplicidad);
    
   //Se comrpueba que la entrada sean solo números.
   regex pattern("^[0-9]+$");
   if (!regex_match(inputGenero, pattern) || !regex_match(inputMultiplicidad, pattern)) {
      cout << "Por favor, introduce valores validos (un unico numero entero sin caracteres ni decimales).\n";
      return 1;
   }
    
   //Cadena a número
   int genero = stoi(inputGenero);
   int multiplicidad = stoi(inputMultiplicidad);
    
   //Condiciones que se deben de cumplir
   if ((genero < 0) || (multiplicidad < 1) || (genero < multiplicidad - 1)) {
      cout << "Por favor, introduce valores validos (genero >= 0, multiplicidad >= 1 y genero >= multiplicidad - 1).\n";
      return 1;
   }
    
   //Condición por defecto: si el género es 0 y la multiplicidad 1
   if (genero == 0 && multiplicidad == 1) {
      cout << "Semigrupos numericos internos:\n<1>\n\nSemigrupos numericos hoja:\n";
      return 0;
   }
    
   cout << "Calculando semigrupos numericos con genero " << genero
      << " y multiplicidad " << multiplicidad << "...\n";
    
   auto inicio = chrono::high_resolution_clock::now();
   encontrarSemigruposYHojas(genero, multiplicidad);
   auto fin = chrono::high_resolution_clock::now();
   auto duracion = chrono::duration_cast<chrono::seconds>(fin - inicio).count();
   
   cout << "\nEl programa tardo " << duracion << " segundos.\n";
   return 0;
}

