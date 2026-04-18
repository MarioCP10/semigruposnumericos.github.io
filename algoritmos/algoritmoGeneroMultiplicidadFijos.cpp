/**
 * @file fixedGenusMultiplicityAlgorithm.cpp
 * @brief Algoritmo para calcular semigrupos numéricos internos y hojas, dado su género y multiplicidad
 * @details Se utiliza OpenMP para paralelizar la búsqueda, cachés divididas (shards) para evitar colisiones 
 * entre hilos, y se captura la señal SIGINT (Ctrl+C) para abortar de forma limpia y mostrar resultados que se lleven
 * hasta ese momento
 */

#include <csignal> 
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
#include <queue>
#include <limits>
#include <deque>
#include <atomic>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads(){ return 1; }
inline int omp_get_thread_num(){ return 0; }
inline void omp_set_num_threads(int){ }
#endif

using namespace std;

/** * @brief Bandera atómica global para indicar si el usuario ha solicitado la interrupción del programa
 */
static std::atomic<bool> interrupted{false};

/**
 * @brief Controlador seguro para la señal SIGINT (Ctrl+C).
 * @details Se cambia el estado de la bandera interrupted a true. Es vital no realizar operaciones de E/S 
 * ni llamadas no reentrantes dentro de esta función
 * @param signum Número de la señal capturada
 */
static void sigint_handler(int) {
  interrupted.store(true);
}

/**
 * @brief Estructura para calcular el hash de un vector de enteros
 * @details Se implementa una función de mezcla para que el vector pueda usarse como clave en un unordered_map
 */
struct VectorHash {
  size_t operator()(const vector<int>& v) const noexcept {
    size_t h = v.size();
    for (int x : v) {
      h ^= std::hash<int>{}(x) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    }
    return h;
  }
};

/** @brief Caché para almacenar resultados del cálculo del género, segmentada por hilo (shard) */
static vector<unordered_map<vector<int>, int, VectorHash>> cacheGeneroShards;

/** @brief Caché para almacenar resultados del cálculo de Frobenius, segmentada por hilo (shard) */
static vector<unordered_map<vector<int>, int, VectorHash>> cacheFrobeniusShards;

/** @brief Límite máximo de elementos por cada caché para evitar el agotamiento de memoria */
const size_t MAX_CACHE_ENTRIES = 200000;

/**
 * @brief Asegura que la caché del hilo especificado no exceda el tamaño máximo
 * @details Si se excede MAX_CACHE_ENTRIES, se limpia la caché para liberar memoria
 * @param shard Identificador del hilo
 */
static inline void ensure_cache_size_for_shard(int shard) {
  if (shard < 0) shard = 0;
  if ((size_t)shard >= cacheGeneroShards.size()) return;
  auto &cg = cacheGeneroShards[shard];
  if (cg.size() > MAX_CACHE_ENTRIES) cg.clear();
  auto &cf = cacheFrobeniusShards[shard];
  if (cf.size() > MAX_CACHE_ENTRIES) cf.clear();
}

/**
 * @brief Calcula el Máximo Común Divisor (MCD) de dos números
 * @param a Primer número
 * @param b Segundo número
 * @return El MCD de a y b
 */
int maxCD(int a, int b) {
  while (b) {
    a %= b;
    swap(a, b);
  }
  return a;
}

/**
 * @brief Comprueba si el MCD de todos los elementos de un conjunto es 1
 * @param S Conjunto de enteros
 * @return true si el MCD global es 1, false en caso contrario
 */
bool mcdEsUno(const vector<int>& S) {
  if (S.empty()) return false;
  int g = S[0];
  for (size_t i = 1; i < S.size() && g != 1; ++i) g = std::gcd(g, S[i]);
  return g == 1;
}

/**
 * @brief Calcula el conjunto de Apéry del semigrupo generado por S con respecto a su elemento mínimo
 * @details Utiliza el algoritmo de Dijkstra sobre las clases de residuos módulo m (donde m = min(S)).
 * Se comprueba si ha habido alguna interrupción
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return Vector w donde w[i] es el elemento más pequeño del semigrupo congruente con i módulo m
 */
static vector<int> apery(const vector<int>& S) {
  int m = *min_element(S.begin(), S.end());
  const int INF = numeric_limits<int>::max() / 4;
  vector<int> w(m, INF);

  priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;
  w[0] = 0;
  pq.emplace(0, 0);

  while (!pq.empty()) {
    if (interrupted.load()) break;
    auto top = pq.top(); pq.pop();
    int cost = top.first;
    int r = top.second;
    if (cost != w[r]) continue;
    for (int a : S) {
      int nr = (r + a) % m;
      int ncost = cost + a;
      if (ncost < w[nr]) {
        w[nr] = ncost;
        pq.emplace(ncost, nr);
      }
    }
  }
  return w;
}

/**
 * @brief Calcula el número de Frobenius de un semigrupo utilizando la caché del hilo en concreto
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @param shard Identificador del hilo actual
 * @return El número de Frobenius, o -1 si el semigrupo contiene todos los enteros positivos
 */
int calculaFrobenius_shard(const vector<int>& S, int shard) {
  vector<int> key = S;
  if (!is_sorted(key.begin(), key.end())) sort(key.begin(), key.end());

  auto &cacheF = cacheFrobeniusShards[shard];
  auto it = cacheF.find(key);
  if (it != cacheF.end()) return it->second;

  auto w = apery(key);
  int m = *min_element(key.begin(), key.end());
  int F = *max_element(w.begin(), w.end()) - m;

  if (F < 0) {
    cacheF[key] = -1;
    return -1;
  }

  cacheF[key] = F;
  return F;
}

/**
 * @brief Wrapper para calcular el número de Frobenius obteniendo automáticamente el hilo actual
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return El número de Frobenius
 */
int calculaFrobenius(const vector<int>& S) {
  int shard = omp_get_thread_num();
  return calculaFrobenius_shard(S, shard);
}

/**
 * @brief Calcula el conductor del semigrupo, siendo este Frobenius + 1.
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return El conductor. Devuelve 0 si Frobenius es -1.
 */
int calculaConductor(const vector<int>& S) {
  int F = calculaFrobenius(S);
  if (F == -1) return 0;
  return F + 1;
}

/**
 * @brief Calcula el género del semigrupo numérico usando el conjunto de Apéry con acceso a caché local
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @param shard Identificador del hilo actual
 * @return El género del semigrupo
 */
int calculaGenero_shard(const vector<int>& S, int shard) {
  vector<int> key = S;
  if (!is_sorted(key.begin(), key.end())) sort(key.begin(), key.end());

  auto &cacheG = cacheGeneroShards[shard];
  auto it = cacheG.find(key);
  if (it != cacheG.end()) return it->second;

  auto w = apery(key);
  int m = *min_element(key.begin(), key.end());
  int F = *max_element(w.begin(), w.end()) - m;
  if (F < 0) {
    cacheG[key] = 0;
    return 0;
  }
  int c = F + 1;

  long long countInSBelowC = 0;
  for (int r = 0; r < m; ++r) {
    if (w[r] <= c - 1) {
      countInSBelowC += 1 + ((c - 1 - w[r]) / m);
    }
  }
  int genero = (int)(c - countInSBelowC);

  cacheG[key] = genero;
  return genero;
}

/**
 * @brief Wrapper para calcular el género obteniendo automáticamente el hilo actual
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return El género
 */
int calculaGenero(const vector<int>& S) {
  int shard = omp_get_thread_num();
  return calculaGenero_shard(S, shard);
}

/**
 * @brief Determina si el conjunto dado forma la base mínima de Hilbert para el semigrupo numérico
 * @details Un generador es redundante si se puede formar combinando los otros. Se comprueba además la bandera de interrupción
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return true si es una base mínima, false si hay redundancia o si se interrumpió el proceso
 */
bool esMinimalHilbert(const vector<int>& S) {
  if (interrupted.load()) return false;
  if (S.size() == 1) 
    return true;

  int m = *min_element(S.begin(), S.end());
  int conductor = calculaConductor(S);
  int limit = conductor + m;

  vector<char> vis(limit + 1, 0);
  deque<int> q;
  vis[0] = 1;
  q.push_back(0);
  while (!q.empty()) {
    if (interrupted.load()) return false;
    int t = q.front(); q.pop_front();
    for (int a : S) {
      int nt = t + a;
      if (nt <= limit && !vis[nt]) {
        vis[nt] = 1;
        q.push_back(nt);
      }
    }
  }

  int consecutivos = 0, prev = -2;
  for (int i = 0; i <= limit; ++i) {
    if (vis[i]) {
       consecutivos = (i == prev + 1) ? consecutivos + 1 : 1;
       if (consecutivos >= m) break;
        prev = i;
    }
  }
  if (consecutivos < m) return false;

  for (size_t i = 0; i < S.size(); ++i) {
    if (interrupted.load()) return false;
    vector<int> copia; copia.reserve(S.size()-1);
    for (size_t j = 0; j < S.size(); ++j) if (j != i) copia.push_back(S[j]);

    int lim = S[i];
    vector<char> vis2(lim + 1, 0);
    deque<int> q2;
    vis2[0] = 1; q2.push_back(0);
    while (!q2.empty()) {
      if (interrupted.load()) return false;
      int t = q2.front(); q2.pop_front();
      for (int a : copia) {
        int nt = t + a;
        if (nt <= lim && !vis2[nt]) {
          vis2[nt] = 1;
          q2.push_back(nt);
        }
      }
    }
    if (vis2[S[i]]) return false;
  }

  for (int s : S) if (!vis[s]) return false;

  return true;
}

/**
 * @brief Comprueba si el semigrupo es una hoja en el árbol de semigrupos para un género fijo
 * @details Una hoja se define como un semigrupo cuyo número de Frobenius es mayor que todos sus generadores mínimos
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @param genero Género esperado
 * @return true si cumple la condición de hoja, false si no o si fue interrumpido
 */
bool esHoja(const vector<int>& S, int genero) {
  if (interrupted.load()) return false;
  if (calculaGenero(S) != genero)
    return false;
  int f = calculaFrobenius(S);
  if (f == -1)
    return false;
  int max_elem = *max_element(S.begin(), S.end());
  return f > max_elem;
}

/**
 * @brief Función recursiva auxiliar para generar combinaciones posibles de generadores
 * @details Realiza una poda del árbol de búsqueda basándose en si es matemáticamente posible 
 * que la rama actual llegue a tener un MCD igual a 1. Cada hilo escribe sus descubrimientos en su buffer local
 * @param numeros Lista de candidatos enteros válidos
 * @param tamano Cantidad de elementos de la combinación actual a buscar
 * @param start Índice de inicio para la combinación actual
 * @param comb Estado actual de la combinación, siendo un acumulador recursivo
 * @param gcd_parcial MCD actual de la combinación temporal
 * @param genero El género fijo
 * @param multiplicidad La multiplicidad requerida
 * @param local_buffer Vector de strings privado del hilo para evitar contención de E/S de la consola
 * @param contador_internos Referencia a contador atómico de semigrupos numéricos internos encontrados
 * @param contador_hojas Referencia a contador atómico de semigrupos numéricos hoja encontrados
 */
void genera_combinaciones_rec(const vector<int>& numeros,
                              int tamano,
                              int start,
                              vector<int>& comb,
                              int gcd_parcial,
                              int genero,
                              int multiplicidad,
                              vector<string>& local_buffer,
                              atomic<size_t>& contador_internos,
                              atomic<size_t>& contador_hojas) {
  if (interrupted.load()) return; // comprobación temprana

  int n = (int)numeros.size();
  int depth = (int)comb.size();
  if (depth == tamano) {
    if (interrupted.load()) return;

    vector<int> S = comb; // Creciente por construcción

    if (*min_element(S.begin(), S.end()) != multiplicidad)
      return;

    if (!mcdEsUno(S)) return;

    int tid = omp_get_thread_num();
    ensure_cache_size_for_shard(tid);

    if (calculaGenero_shard(S, tid) == genero) {
      if (esMinimalHilbert(S)) {
        string out = "<";
        for (size_t i = 0; i < S.size(); ++i) {
          out += to_string(S[i]);
          if (i + 1 < S.size()) out += ",";
        }
        out += ">\n";

        // Almacenamos en el buffer local del hilo tanto internos como hojas
        local_buffer.push_back(out);

        if (esHoja(S, genero)) ++contador_hojas;
        else ++contador_internos;
      }
    }
    return;
  }

  for (int i = start; i <= n - (tamano - depth); ++i) {
    if (interrupted.load()) return;

    int val = numeros[i];
    int new_gcd = (depth == 0) ? val : std::gcd(gcd_parcial, val);

    if (new_gcd != 1) {
      bool possible = false;
      for (int cand = i + 1; cand < n; ++cand) {
          if (std::gcd(new_gcd, numeros[cand]) == 1) { possible = true; break; }
      }
      if (!possible) continue;
    }

    comb.push_back(val);
    genera_combinaciones_rec(numeros, tamano, i + 1, comb, new_gcd, genero, multiplicidad,
                             local_buffer, contador_internos, contador_hojas);
    comb.pop_back();

    if (interrupted.load()) return;
  }
}

/**
 * @brief Bucle principal que organiza la búsqueda y clasificación de semigrupos numéricos en paralelo
 * @details Se prereserva memoria, se configura el entorno paralelo, disparando la recursión con sharding 
 * y centralizando y formateando la salida de resultados y conteos
 * @param genero El género fijo proporcionado por el usuario
 * @param multiplicidad La multiplicidad fija proporcionada por el usuario
 */
void encontrarSemigruposYHojas(int genero, int multiplicidad) {
  std::atomic<size_t> contador_internos{0};
  std::atomic<size_t> contador_hojas{0};
  
  // Espacio de búsqueda
  int limite = (genero * 3) - (genero - 1);
  if (limite < multiplicidad) limite = multiplicidad;
  vector<int> numeros;
  for (int i = multiplicidad; i <= limite; ++i)
    numeros.push_back(i);

  cout << "\nSemigrupos numericos internos (m=" << multiplicidad << ", g=" << genero << "):\n";

  int n = (int)numeros.size();

  int desired_threads = 8;
  #ifdef _OPENMP
    omp_set_num_threads(desired_threads);
  #endif

  int max_threads = omp_get_max_threads();
  if (max_threads <= 0) max_threads = 1;
  cout << "Solicitados " << desired_threads << " hebras. OpenMP reporta max_threads = " << max_threads << ".\n";

  cacheGeneroShards.assign(max_threads, unordered_map<vector<int>, int, VectorHash>());
  cacheFrobeniusShards.assign(max_threads, unordered_map<vector<int>, int, VectorHash>());

  for (int tamano = 2; tamano <= genero; ++tamano) {
    if (interrupted.load()) break;

    int i0_max = n - tamano;
    if (i0_max < 0) break;

    vector<vector<string>> buffers_por_hilo(max_threads);

    #pragma omp parallel for schedule(dynamic)
    for (int i0 = 0; i0 <= i0_max; ++i0) {
      if (interrupted.load()) continue;
      int tid = omp_get_thread_num();
      if (tid < 0) tid = 0;
      tid = tid % max_threads; // Nos aseguramos que se obtengan los semigrupos que cumplan las condiciones especificadas

      auto &local_buffer = buffers_por_hilo[tid];

      vector<int> comb;
      comb.reserve(tamano);
      comb.push_back(numeros[i0]);
      int gcd_parcial = numeros[i0];

      genera_combinaciones_rec(numeros, tamano, i0 + 1, comb, gcd_parcial, genero, multiplicidad,
                                  local_buffer, contador_internos, contador_hojas);
    }

    // Volcamos buffers por hilo en orden (esto imprimirá lo ya encontrado)
    for (int t = 0; t < max_threads; ++t) {
      for (auto &s : buffers_por_hilo[t]) {
         cout << s;
      }
    }
    // Si recibimos interrupción, no continuaremos con siguientes tamaños
    if (interrupted.load()) break;
   }

   // Para reducir el espacio de búsqueda, añadimos el semigrupo que siempre se contiene para el genero dado y la correspondiente multiplicidad
   if (!interrupted.load() && multiplicidad == genero + 1) {
     vector<int> semigrupoExtra;
     for (int i = genero + 1; i <= 2 * genero + 1; ++i)
       semigrupoExtra.push_back(i);
     cout << "<";
     for (size_t i = 0; i < semigrupoExtra.size(); ++i)
       cout << semigrupoExtra[i] << (i + 1 < semigrupoExtra.size() ? "," : "");
     cout << ">\n";
     ++contador_internos;
   }

   // Imprimimos semigrupos hojas (fase secuencial para reproducibilidad)
   if (!interrupted.load()) {
     for (int tamano = 2; tamano <= genero; ++tamano) {
         if (interrupted.load()) break;
           vector<int> comb;
           function<void(int,int)> backtrack;
           backtrack = [&](int start, int depth){
             if (interrupted.load()) return;
             if ((int)comb.size() == tamano) {
                vector<int> S = comb;
                if (*min_element(S.begin(), S.end()) != multiplicidad) return;
                if (!mcdEsUno(S)) return;
                int tid = 0;
                ensure_cache_size_for_shard(tid);
                if (calculaGenero_shard(S, tid) == genero) {
                  if (esMinimalHilbert(S) && esHoja(S, genero)) {
                     cout << "<";
                     for (size_t i = 0; i < S.size(); ++i)
                        cout << S[i] << (i + 1 < S.size() ? "," : "");
                     cout << ">\n";
                  }
                }
                return;
             }
             int n_local = (int)numeros.size();
             int need = tamano - (int)comb.size();
             for (int i = start; i <= n_local - need; ++i) {
               if (interrupted.load()) return;
                 comb.push_back(numeros[i]);
                 backtrack(i + 1, depth + 1);
                 comb.pop_back();
             }
           };
      backtrack(0, 0);
     }
   } 
   else {
      cout << "\n*** Interrupción recibida: se muestran resultados parciales hasta el momento. ***\n";
   }

   // Conteo final (fase secuencial para consistencia)
  size_t count_internos_final = 0;
  size_t count_hojas_final = 0;

  if (!interrupted.load()) {
    for (int tamano = 2; tamano <= genero; ++tamano) {
      if (interrupted.load()) break;
        vector<int> comb;
        function<void(int)> backtrack_count;
        backtrack_count = [&](int start) {
          if (interrupted.load()) return;
          if ((int)comb.size() == tamano) {
             vector<int> S = comb;
              if (*min_element(S.begin(), S.end()) != multiplicidad) return;
              if (!mcdEsUno(S)) return;
              int tid = 0;
              ensure_cache_size_for_shard(tid);
              if (calculaGenero_shard(S, tid) == genero && esMinimalHilbert(S)) {
                  if (esHoja(S, genero)) ++count_hojas_final;
                  else ++count_internos_final;
              }
              return;
          }
          int n_local = (int)numeros.size();
          int need = tamano - (int)comb.size();
          for (int i = start; i <= n_local - need; ++i) {
             if (interrupted.load()) return;
               comb.push_back(numeros[i]);
               backtrack_count(i + 1);
               comb.pop_back();
          }
        };
        backtrack_count(0);
    }

    if (multiplicidad == genero + 1) ++count_internos_final;
   } 
   else {
     cout << "Conteos finales no completos debido a la interrupción; mostrando conteos parciales si existen.\n";
     // Podemos usar los contadores atómicos incrementados durante la fase paralela
     count_internos_final = contador_internos.load();
     count_hojas_final = contador_hojas.load();
   }

  cout << "\nSemigrupos numericos hoja (m=" << multiplicidad << ", g=" << genero << "):\n";
  cout << "\nCantidad de internos: " << count_internos_final << "\n\n";
  cout << "Cantidad de hojas: " << count_hojas_final << "\n";

  if (count_internos_final > count_hojas_final)
    cout << "\nComparacion:\nHay mas semigrupos numericos internos que hojas.\n";
  else if (count_internos_final < count_hojas_final)
    cout << "\nComparacion:\nHay mas semigrupos numericos hojas que internos.\n";
  else
    cout << "\nComparacion:\nHay igual cantidad de semigrupos numericos internos y hojas.\n";
}

/**
 * @brief Función principal del programa
 * @details Se configura la captura de la señal de interrupción, se solicita los inputs básicos al usuario, 
 * se valida los datos de entrada y se ejecuta la búsqueda cronometrando el tiempo de cómputo
 * @return 0 en caso de éxito, 1 ante entradas inválidas o excepciones capturadas.
 */
int main() {
  try {
    // Registramos handler de SIGINT
    std::signal(SIGINT, sigint_handler);

    string inputGenero, inputMultiplicidad;

    cout << "Introduce el genero: ";
    getline(cin, inputGenero);

    cout << "Introduce la multiplicidad: ";
    getline(cin, inputMultiplicidad);

    regex pattern("^[0-9]+$");
    if (!regex_match(inputGenero, pattern) || !regex_match(inputMultiplicidad, pattern)) {
       cout << "Por favor, introduce valores validos (un unico numero entero sin caracteres ni decimales).\n";
       return 1;
    }

    int genero = stoi(inputGenero);
    int multiplicidad = stoi(inputMultiplicidad);

    if ((genero < 0) || (multiplicidad < 1) || (genero < multiplicidad - 1)) {
     cout << "Por favor, introduce valores validos (genero >= 0, multiplicidad >= 1 y genero >= multiplicidad - 1).\n";
     return 1;
    }

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
   } 
   catch (const std::exception& e) {
      cerr << "Excepción atrapada: " << e.what() << endl;
      return 1;
   } 
   catch (...) {
      cerr << "Excepción desconocida atrapada.\n";
      return 1;
   }

  return 0;
}

