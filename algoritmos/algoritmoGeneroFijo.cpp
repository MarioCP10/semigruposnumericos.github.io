/**
 * @file fixedGenusAlgorithm.cpp
 * @brief Programa para calcular la cantidad de semigrupos numéricos existentes a partir de un género dado. Estos semigrupos se clasifican como internos u hojas
 * @details Se utiliza OpenMP para paralelización, algoritmo de Dijkstra para el conjunto de Apéry y un sistema de cachés particionadas por hilo (sharding) para optimizar el rendimiento del programa
 */

#include <bits/stdc++.h>
#include <omp.h>
#include <atomic>
using namespace std;

/**
 * @brief Estructura para calcular el hash de un vector de enteros
 * @details Se utiliza como función de hash para almacenar vectores como claves en un unordered_map
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

/** @brief Número máximo de entradas permitidas en la caché antes de limpiarla */
const size_t MAX_CACHE_ENTRIES = 200000;

// Caches por shard (dos modos: mask o vector)

/** @brief Caché particionada por hilo (shard) para el género usando vector como clave */
static vector<unordered_map<vector<int>, int, VectorHash>> cacheGeneroShards_vec;
/** @brief Caché particionada por hilo (shard) para el número de Frobenius usando vector como clave */
static vector<unordered_map<vector<int>, int, VectorHash>> cacheFrobeniusShards_vec;

/** @brief Caché particionada por hilo (shard) para el género usando máscara de bits (uint64_t) como clave */
static vector<unordered_map<uint64_t, int>> cacheGeneroShards_mask;
/** @brief Caché particionada por hilo (shard) para el número de Frobenius usando máscara de bits (uint64_t) como clave */
static vector<unordered_map<uint64_t, int>> cacheFrobeniusShards_mask;

// Configuración global para caches y array de numeros

/** @brief Indica si se debe usar la caché basada en máscaras de bits, siendo más rápida para números pequeños */
static bool use_mask_cache = false;
/** @brief Base global para la máscara de bits, donde el primer número corresponde al bit 0 en modo mask */
static int global_base = 2;
/** @brief Lista global de números a combinar */
static vector<int> global_numeros;
/** @brief Vector que almacena el máximo común divisor de los sufijos para podar recursiones de forma anticipada */
static vector<int> global_suffix_gcd;
/** @brief Tamaño global de la lista de números */
static int global_n = 0;

/**
 * @brief Asegura que la caché del hilo especificado no supere el límite máximo
 * @details Si supera MAX_CACHE_ENTRIES, se limpia para evitar un consumo excesivo de memoria
 * @param shard El identificador del hilo (shard) actual
 */
static inline void ensure_cache_size_for_shard(int shard) {
  if (use_mask_cache) {
    auto &cg = cacheGeneroShards_mask[shard];
    if (cg.size() > MAX_CACHE_ENTRIES) cg.clear();
    auto &cf = cacheFrobeniusShards_mask[shard];
    if (cf.size() > MAX_CACHE_ENTRIES) cf.clear();
  } 
  else {
    auto &cg = cacheGeneroShards_vec[shard];
    if (cg.size() > MAX_CACHE_ENTRIES) cg.clear();
    auto &cf = cacheFrobeniusShards_vec[shard];
    if (cf.size() > MAX_CACHE_ENTRIES) cf.clear();
  }
}

/**
 * @brief Calcula el Máximo Común Divisor de dos números
 * @param a Primer número
 * @param b Segundo número
 * @return El Máximo Común Divisor de a y b
 */
int maxCD(int a, int b) {
  while (b) {
    a %= b;
    std::swap(a, b);
  }
  return a;
}

/**
 * @brief Comprueba si el Máximo Común Divisor de un conjunto de números es 1
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return true si el MCD es 1, false en caso contrario o si está vacío
 */
bool mcdEsUno(const vector<int>& S) {
  if (S.empty()) return false;
  int g = S[0];
  for (size_t i = 1; i < S.size() && g != 1; ++i) g = std::gcd(g, S[i]);
  return g == 1;
}

/**
 * @brief Calcula el conjunto de Apéry de un semigrupo con respecto a su multiplicidad
 * @details Se utiliza el algoritmo de Dijkstra sobre los residuos módulo m
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return Un vector que representa el conjunto de Apéry w(S)
 */
static vector<int> apery(const vector<int>& S) {
  int m = *min_element(S.begin(), S.end());
  const int INF = numeric_limits<int>::max() / 4;
  vector<int> w(m, INF);

  // Dijkstra en grafos de m nodos
  priority_queue<pair<int,int>, vector<pair<int,int>>,
                 greater<pair<int,int>>> pq;
  w[0] = 0;
  pq.emplace(0, 0);

  while (!pq.empty()) {
    auto [cost, r] = pq.top(); pq.pop();
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
 * @brief Crea una clave basada en máscara de bits a partir de los generadores
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return Clave de 64 bits representando el conjunto
 */
static inline uint64_t make_mask_key(const vector<int>& S) {
  uint64_t key = 0;
  for (int v : S) key |= (1ULL << (v - global_base));
  return key;
}

/**
 * @brief Calcula el número de Frobenius accediendo a la caché del hilo específico
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @param shard ID del hilo de OpenMP
 * @return El número de Frobenius del semigrupo numérico
 */
int calculaFrobenius_shard(const vector<int>& S, int shard) {
  if (use_mask_cache) {
    uint64_t key = make_mask_key(S);
    auto &cacheF = cacheFrobeniusShards_mask[shard];
    auto it = cacheF.find(key);
    if (it != cacheF.end()) return it->second;

    auto w = apery(S);
    int m = *min_element(S.begin(), S.end());
    int F = *max_element(w.begin(), w.end()) - m;

    cacheF.emplace(key, F);
    return F;
  } 
  else {
    auto &cacheF = cacheFrobeniusShards_vec[shard];
    auto it = cacheF.find(S);
    if (it != cacheF.end()) return it->second;

    auto w = apery(S);
    int m = *min_element(S.begin(), S.end());
    int F = *max_element(w.begin(), w.end()) - m;

    cacheF.emplace(S, F);
    return F;
  }
}

/**
 * @brief Wrapper para calcular el número de Frobenius obteniendo automáticamente el hilo actual
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return El número de Frobenius del semigrupo
 */
int calculaFrobenius(const vector<int>& S) {
  int shard = omp_get_thread_num();
  return calculaFrobenius_shard(S, shard);
}

/**
 * @brief Calcula el conductor del semigrupo numérico
 * @details El conductor es siempre el número de Frobenius + 1
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return El conductor del semigrupo numérico
 */
int calculaConductor(const vector<int>& S) {
  int F = calculaFrobenius(S);
  return F + 1;
}

/**
 * @brief Calcula el género de un semigrupo y guarda en caché sus resultados
 * @details Se aprovecha el cálculo para almacenar también el número de Frobenius en caché
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @param shard ID del hilo actual
 * @return El género del semigrupo numérico
 */
int calculaGenero_shard(const vector<int>& S, int shard) {
  if (use_mask_cache) {
    uint64_t key = make_mask_key(S);
    auto &cacheG = cacheGeneroShards_mask[shard];
    auto it = cacheG.find(key);
    if (it != cacheG.end()) return it->second;

    auto w = apery(S);
    int m = *min_element(S.begin(), S.end());
    int F = *max_element(w.begin(), w.end()) - m;
    int c = F + 1;

    long long countInSBelowC = 0;
    for (int r = 0; r < m; ++r) {
      if (w[r] <= c - 1) {
        countInSBelowC += 1 + ((c - 1 - w[r]) / m);
      }
    }
    int genero = (int)(c - countInSBelowC);

    cacheG.emplace(key, genero);
    // Cacheamos F en la tabla de Frobenius para evitar recomputos posteriores
    cacheFrobeniusShards_mask[shard].emplace(key, F);
    return genero;
  } 
  else {
    auto &cacheG = cacheGeneroShards_vec[shard];
    auto it = cacheG.find(S);
    if (it != cacheG.end()) return it->second;

    auto w = apery(S);
    int m = *min_element(S.begin(), S.end());
    int F = *max_element(w.begin(), w.end()) - m;
    int c = F + 1;

    long long countInSBelowC = 0;
    for (int r = 0; r < m; ++r) {
      if (w[r] <= c - 1) {
        countInSBelowC += 1 + ((c - 1 - w[r]) / m);
      }
    }
    int genero = (int)(c - countInSBelowC);

    cacheG.emplace(S, genero);
    // Cacheamos F en la tabla de Frobenius para evitar recomputos posteriores
    cacheFrobeniusShards_vec[shard].emplace(S, F);
    return genero;
  }
}

/**
 * @brief Wrapper para calcular el género obteniendo automáticamente el hilo actual
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return El género del semigrupo
 */
int calculaGenero(const vector<int>& S) {
    int shard = omp_get_thread_num();
    return calculaGenero_shard(S, shard);
}

/**
 * @brief Comprueba si un elemento pertenece al semigrupo numérico usando el conjunto de Apéry
 * @param w El conjunto de Apéry calculado
 * @param m La multiplicidad
 * @param x El valor a comprobar
 * @return true si pertenece, false en caso contrario
 */
inline bool enSemigrupo_porApery(const vector<int>& w, int m, int x) {
  return x >= w[x % m];
}

/**
 * @brief Comprueba si el conjunto de generadores es minimal a partir de la base de Hilbert
 * @details Se comprueba que ningún generador pueda ser expresado como combinación lineal del resto
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @return true si es minimal, false si existen generadores redundantes
 */
bool esMinimalHilbert(const vector<int>& S) {
  if (S.size() == 1) return true;

  int m = *min_element(S.begin(), S.end());
  int conductor = calculaConductor(S); // Estará en cache si ya se calculó el genero
  int limit = conductor + m;

  // BFS para cierre hasta limit
  vector<char> vis(limit + 1, 0);
  deque<int> q;
  vis[0] = 1;
  q.push_back(0);
  while (!q.empty()) {
    int t = q.front(); q.pop_front();
    for (int a : S) {
      int nt = t + a;
      if (nt <= limit && !vis[nt]) {
        vis[nt] = 1;
        q.push_back(nt);
      }
    }
  }

  // Bloque consecutivo de longitud m
  int consecutivos = 0, prev = -2;
  for (int i = 0; i <= limit; ++i) {
    if (vis[i]) {
      consecutivos = (i == prev + 1) ? consecutivos + 1 : 1;
      if (consecutivos >= m) break;
      prev = i;
    }
  }
  
  if (consecutivos < m) return false;

  // Minimalidad: cada generador no debe representarse con los demás
  // Observación: comprobar si a_i ∈ <S\{a_i}>. Mantenemos el método BFS sobre [0..a_i]
  for (size_t i = 0; i < S.size(); ++i) {
    vector<int> copia; copia.reserve(S.size()-1);
    for (size_t j = 0; j < S.size(); ++j) if (j != i) copia.push_back(S[j]);

    int lim = S[i];
    vector<char> vis2(lim + 1, 0);
    deque<int> q2;
    vis2[0] = 1; 
    q2.push_back(0);
    while (!q2.empty()) {
      int t = q2.front(); 
      q2.pop_front();
      for (int a : copia) {
        int nt = t + a;
        if (nt <= lim && !vis2[nt]) {
          vis2[nt] = 1;
          q2.push_back(nt);
        }
      }
    }
    if (vis2[S[i]]) return false; // No es minimal
  }

  for (int s : S) if (!vis[s]) return false;

  return true;
}

/**
 * @brief Determina si un semigrupo es una hoja en el árbol de semigrupos numéricos
 * @details Un semigrupo es hoja si su número de Frobenius es mayor que todos sus generadores mínimos
 * @param S Vector de enteros que representa los generadores del semigrupo numérico
 * @param genero Género fijo a comprobar
 * @return true si es hoja, false si es nodo interno o no coincide el género
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
 * @brief Genera combinaciones recursivamente, podando ramas inviables
 * @details Se construyen posibles conjuntos de generadores y se comprueba su género y minimalidad si el MCD parcial permite que llegue a 1
 * @param numeros Vector de candidatos
 * @param tamano Tamaño de la combinación buscada
 * @param start Índice inicial para iterar, evitando duplicados
 * @param comb Vector temporal que almacena la combinación actual
 * @param gcd_parcial El MCD de la combinación parcial
 * @param genero Género fijo de los semigrupos
 * @param local_buffer Búfer local del hilo para acumular resultados que se imprimen
 * @param contador_internos Referencia atómica al contador global de nodos internos
 * @param contador_hojas Referencia atómica al contador global de nodos hojas
 */
void genera_combinaciones_rec(const vector<int>& numeros, int tamano, int start, vector<int>& comb, int gcd_parcial, int genero, vector<string>& local_buffer, atomic<size_t>& contador_internos, atomic<size_t>& contador_hojas) {
  int n = (int)numeros.size();
  int depth = (int)comb.size();
  if (depth == tamano) {
    // comb ya está ordenado y sin duplicados por la generación
    if (!mcdEsUno(comb)) return;

    int tid = omp_get_thread_num();
    ensure_cache_size_for_shard(tid);

    if (calculaGenero_shard(comb, tid) == genero) {
      if (esMinimalHilbert(comb)) {
        string out = "<";
        for (size_t i = 0; i < comb.size(); ++i) {
          out += to_string(comb[i]);
          if (i + 1 < comb.size()) out += ",";
        }
        out += ">\n";
        local_buffer.push_back(std::move(out));

        if (esHoja(comb, genero)) ++contador_hojas;
        else ++contador_internos;
      }
    }
    return;
  }

  // Poda por límites. Si no hay suficientes elementos restantes
  for (int i = start; i <= n - (tamano - depth); ++i) {
    int val = numeros[i];
    int new_gcd = (depth == 0) ? val : std::gcd(gcd_parcial, val);

    // Poda por Máximo Común Divisor usando suffix_gcd. Si incluso con todos los elementos restantes no se puede alcanzar gcd 1, lo saltamos
    if (new_gcd != 1) {
      if (i + 1 >= n) continue;
      int suffix = global_suffix_gcd[i + 1];
      if (std::gcd(new_gcd, suffix) != 1) continue;
    }

    comb.push_back(val);
    genera_combinaciones_rec(numeros, tamano, i + 1, comb, new_gcd, genero, local_buffer, contador_internos, contador_hojas);
    comb.pop_back();
  }
}

/**
 * @brief Método de búsqueda paralela de semigrupos dado un género
 * @details Se inicializa cachés, buffers, se genera el espacio de búsqueda y se inicia la generación recursiva con OpenMP.
 * @param genero El género de los semigrupos que se quieren encontrar
 */
void encontrarSemigruposYHojas(int genero) {
  std::atomic<size_t> contador_internos{0};
  std::atomic<size_t> contador_hojas{0};

  // Espacio de búsqueda
  int limite = (genero * 3) - (genero - 1);
  if (limite < 3) limite = 3;
  vector<int> numeros(limite - 1);
  iota(numeros.begin(), numeros.end(), 2);

  cout << "\nSemigrupos numericos internos y hojas:\n";

  int n = (int)numeros.size();

  // Preparamos caches por hilo y decidimos modo mask o vector
  int max_threads = omp_get_max_threads();

  global_numeros = numeros;
  global_n = n;
  global_base = 2; // Los primeros numeros empiezan en 2

  use_mask_cache = (n <= 64); // Modo mask si cabe en 64 bits
  if (use_mask_cache) {
    cacheGeneroShards_mask.assign(max_threads, unordered_map<uint64_t,int>());
    cacheFrobeniusShards_mask.assign(max_threads, unordered_map<uint64_t,int>());
    for (int t = 0; t < max_threads; ++t) {
      cacheGeneroShards_mask[t].reserve(1024);
      cacheFrobeniusShards_mask[t].reserve(1024);
    }
  } 
  else {
    cacheGeneroShards_vec.assign(max_threads, unordered_map<vector<int>,int,VectorHash>());
    cacheFrobeniusShards_vec.assign(max_threads, unordered_map<vector<int>,int,VectorHash>());
    for (int t = 0; t < max_threads; ++t) {
      cacheGeneroShards_vec[t].reserve(1024);
      cacheFrobeniusShards_vec[t].reserve(1024);
    }
  }

  // Construir suffix_gcd global para poda por MCD
  global_suffix_gcd.assign(n, 0);
  if (n > 0) {
    global_suffix_gcd[n-1] = numeros[n-1];
    for (int i = n - 2; i >= 0; --i) global_suffix_gcd[i] = std::gcd(numeros[i], global_suffix_gcd[i+1]);
  }

  for (int tamano = 2; tamano <= genero; ++tamano) {
    int i0_max = n - tamano;
    if (i0_max < 0) break;

    vector<vector<string>> buffers_por_hilo(max_threads);

    #pragma omp parallel for schedule(dynamic)
    for (int i0 = 0; i0 <= i0_max; ++i0) {
      int tid = omp_get_thread_num();
      auto &local_buffer = buffers_por_hilo[tid];

      vector<int> comb;
      comb.reserve(tamano);

      // Empezamos la recursión con el primer elemento fijado en i0
      comb.push_back(numeros[i0]);
      int gcd_parcial = numeros[i0];

      // Llamada recursiva que continúa desde i0+1
      genera_combinaciones_rec(numeros, tamano, i0 + 1, comb, gcd_parcial, genero, local_buffer, contador_internos, contador_hojas);
    }

    // Volcamos buffers por hilo en orden
    for (int t = 0; t < max_threads; ++t) {
      for (auto &s : buffers_por_hilo[t]) cout << s;
    }
  }

  // Semigrupo que siempre cumple la condición de género fijo, por lo que lo obtenemos sin aumentar el espacio de búsqueda
  vector<int> semigrupoExtra;
  for (int i = genero + 1; i <= 2 * genero + 1; ++i)
    semigrupoExtra.push_back(i);
  cout << "<";
  for (size_t i = 0; i < semigrupoExtra.size(); ++i)
    cout << semigrupoExtra[i] << (i + 1 < semigrupoExtra.size() ? "," : "");
  cout << ">\n";
  ++contador_internos;

  cout << "\nCantidad de internos: " << contador_internos.load() << "\n\n";
  cout << "Cantidad de hojas: " << contador_hojas.load() << "\n";

  if (contador_internos.load() > contador_hojas.load())
    cout << "\nComparacion:\nHay mas semigrupos numericos internos que hojas.\n";
  else if (contador_internos.load() < contador_hojas.load())
    cout << "\nComparacion:\nHay mas semigrupos numericos hojas que internos.\n";
  else
    cout << "\nComparacion:\nHay igual cantidad de semigrupos numericos internos y hojas.\n";
}

/**
 * @brief Función principal del programa
 * @details Se le pide al usuario el género, se valida la entrada y se ejecuta el cálculo de tiempos
 * @return 0 si finalizó con éxito, 1 en caso de error o entrada inválida
 */
int main() {
  try {
    int genero;
    cout << "Introduce el genero: ";

    string input;
    getline(cin, input);

    regex patron("^[0-9]+$");
    if (!regex_match(input, patron)) {
      cout << "Entrada no valida. Debes introducir un unico numero entero (0 o mayor) sin caracteres." << endl;
      return 1;
    }

    try {
      genero = stoi(input);
    } 
    catch (...) {
      cout << "No se pudo convertir la entrada. Asegúrate de ingresar un número válido." << endl;
      return 1;
    }

    if (genero == 0) {
      cout << "Semigrupos numericos internos:\n";
      cout << "<1>\n";
      cout << "Cantidad de internos: 1\n";
      cout << "\nSemigrupos numericos hoja:\n";
      cout << "Cantidad de hojas: 0\n";
      return 0;
    }

    cout << "Calculando semigrupos numericos internos y hojas...\n";

    auto inicio = chrono::high_resolution_clock::now();
    encontrarSemigruposYHojas(genero);
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

