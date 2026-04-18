/**
 * @file fixedFrobeniusAlgorithm.cpp
 * @brief Generación de semigrupos numéricos dado un número de Frobenius fijo
 * @details Se utiliza representación de estados mediante palabras de 64 bits, una tabla de dispersión personalizada
 * con índices empaquetados a 24 bits para reducir el consumo de RAM y OpenMP para procesamiento paralelo basado en sharding
 */

#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using uint16 = uint16_t;
using uint32 = uint32_t;
using uint64 = uint64_t;
using int64 = long long;

// Configuración

/** @brief Factor de carga máximo para la tabla hash antes de hacer un rehash */
static const double MAX_LOAD_FACTOR = 0.62;

/** @brief Estimación de la cantidad de estados principales, el cual se ajusta para prereservar memoria y evitar reallocs */
size_t expected_states_main = 0;

/** @brief Flag para activar una compactación agresiva si la memoria está muy ajustada */
static const bool AGGRESSIVE_COMPACT = true;

// DP bitset por palabras

/**
 * @brief Estructura para programación dinámica basada en bits
 * @details Permite representar de forma eficiente grandes conjuntos de booleanos empaquetados en palabras de 64 bits
 */
struct BitDP {
  vector<uint64> W; //Vector que almacena los bits en bloques de 64
  
  /**
   * @brief Ajusta el tamaño de la estructura para acomodar al menos nbits
   * @param nbits Número de bits que se necesitan
   */
  void resize_bits(size_t nbits) {
    size_t nwords = (nbits + 63) >> 6;
    W.assign(nwords, 0ULL);
  }
  
  /**
   * @brief Enciende (pone a 1) el bit en la posición indicada
   * @param i Índice del bit a encender
   */
  inline void setbit(size_t i) { W[i >> 6] |= (1ULL << (i & 63)); }
  
  /**
   * @brief Comprueba si el bit en la posición indicada está encendido
   * @param i Índice del bit a comprobar
   * @return true si el bit es 1, false si es 0
   */
  inline bool testbit(size_t i) const { return (W[i >> 6] >> (i & 63)) & 1ULL; }
};

/**
 * @brief Comprueba si un valor es representable como combinación lineal no negativa de los generadores del semigrupo
 * @param valor El número a comprobar
 * @param generadores Lista de generadores del semigrupo
 * @return true si el valor pertenece al semigrupo, false en caso contrario
 */
bool esRepresentable(int valor, const vector<uint16>& generadores) {
  if (valor < 0) return false;
  if (valor == 0) return true;
  size_t nbits = (size_t)valor + 1;
  BitDP dp; dp.resize_bits(nbits);
  dp.setbit(0);
  for (size_t i = 0; i <= (size_t)valor; ++i) {
    if (!dp.testbit(i)) continue;
    for (uint16 g : generadores) {
      size_t nx = i + (size_t)g;
      if (nx <= (size_t)valor) dp.setbit(nx);
    }
  }
  return dp.testbit((size_t)valor);
}

/**
 * @brief Elimina generadores redundantes para obtener la base de Hilbert mínima
 * @param generadores Lista de generadores del semigrupo
 * @return Un nuevo vector estrictamente con los generadores mínimos
 */
vector<uint16> minimizarGeneradores(const vector<uint16>& generadores) {
  vector<uint16> g = generadores;
  sort(g.begin(), g.end());
  vector<uint16> minimal; minimal.reserve(g.size());
  for (uint16 v : g) {
    if (minimal.empty()) { minimal.push_back(v); continue; }
    if (!esRepresentable(v, minimal)) minimal.push_back(v);
  }
  return minimal;
}

/**
 * @brief Comprueba que un número concreto actúe correctamente como número de Frobenius
 * @details El número de Frobenius F no debe ser representable, pero F+1 sí debe de serlo
 * @param generadores Lista de generadores del semigrupo
 * @param F El número de Frobenius válido
 * @return true si F es el número de Frobenius, false en caso contrario
 */
bool frobeniusValido(const vector<uint16>& generadores, int F) {
  return !esRepresentable(F, generadores) && esRepresentable(F + 1, generadores);
}

/**
 * @brief Genera el semigrupo inicial o la hoja madre dado un número de Frobenius
 * @details Se crea un semigrupo compuesto por todos los enteros en el rango [F+1, 2F+1]
 * @param F El número de Frobenius válido
 * @return Vector con los generadores iniciales
 */
vector<uint16> semigrupoInicial(int F) {
  vector<uint16> g; g.reserve(F+1);
  for (int i = F+1; i <= 2*F+1; ++i) g.push_back((uint16)i);
  return g;
}

/**
 * @brief Convierte un vector de generadores en un array de palabras de 64 bits, siendo la máscara
 * @param g Vector de generadores del semigrupo
 * @param words Puntero al array de destino donde se guardarán los bits
 * @param nwords Cantidad máxima de palabras disponibles en el array destino
 * @param max_bit El bit más alto permitido, que normalmente será 2*F + 1
 */
void genvec_to_words(const vector<uint16>& g, uint64* words, size_t nwords, unsigned max_bit) {
  for (size_t i=0;i<nwords;++i) words[i]=0ULL;
  for (uint16 v : g) {
    if ((unsigned)v <= max_bit) {
      size_t wi = ((size_t)v)>>6;
      if (wi >= nwords) continue;
      unsigned bi = v & 63;
      words[wi] |= (1ULL << bi);
    }
  }
}

/**
 * @brief Convierte una lista cruda de generadores en un array de palabras de 64 bits
 * @param gens Puntero a la lista de generadores
 * @param len Longitud de la lista de generadores donde se guardarán los bits
 * @param words Puntero al array de destino
 * @param nwords Cantidad máxima de palabras disponibles en el array destino
 * @param max_bit El bit más alto permitido
 */
void genlist_to_words(const uint16_t* gens, size_t len, uint64* words, size_t nwords, unsigned max_bit) {
  for (size_t i = 0; i < nwords; ++i) words[i] = 0ULL;
  for (size_t i = 0; i < len; ++i) {
    uint16 v = gens[i];
    if ((unsigned)v <= max_bit) {
      size_t wi = ((size_t)v)>>6;
      if (wi >= nwords) continue;
      unsigned bi = v & 63;
      words[wi] |= (1ULL << bi);
    }
  }
}

/**
 * @brief Convierte un array de palabras de 64 bits de vuelta a un vector de generadores enteros
 * @details Se utiliza intrínsecos del compilador __builtin_ctzll cuando están disponibles para contar ceros finales y encontrar los bits encendidos muy rápido
 * @param words Puntero al array de destino
 * @param nwords Cantidad máxima de palabras disponibles en el array destino
 * @return Un vector ordenado con los generadores obtenidos
 */
vector<uint16> words_to_genvec(const uint64* words, size_t nwords) {
  vector<uint16> out; out.reserve(16);
  for (size_t wi = 0; wi < nwords; ++wi) {
    uint64 w = words[wi];
    while (w) {
      #if defined(__GNUG__) || defined(__clang__)
        int tz = __builtin_ctzll(w);
        unsigned idxv = (unsigned)(wi*64 + tz);
        out.push_back((uint16)idxv);
        w &= (w-1);
      #else
        for (int b = 0; b < 64; ++b) if ((w>>b)&1ULL) out.push_back((uint16)(wi*64 + b));
        break;
      #endif
    }
  }
  sort(out.begin(), out.end());
  return out;
}

/**
 * @brief Estructura de tabla compacta que empaqueta índices en 3 bytes en formato little-endian
 * @details Se guarda el índice en 3 bytes en vez de uint32_t, que ocupa 4 bytes. Esto ahorra 25% de memoria en la tabla de buckets
 * El valor 0 indica que la posición está vacía. Un valor N almacenado corresponde al índice N-1
 */
struct Packed24Table {
  vector<uint8_t> data; // Array lineal subyacente. Tamaño = table_size * 3
  size_t table_size;    // Cantidad de buckets lógicos en la tabla. Debe ser potencia de 2
  size_t mask;          // Máscara de bits (table_size - 1) para operaciones rápidas de módulo
  
  Packed24Table(): table_size(0), mask(0) {}
  
  /**
   * @brief Inicializa la tabla reservando espacio para al menos ts elementos
   * @param ts Tamaño deseado de la tabla, el cual se redondeará a la siguiente potencia de 2
   */
  void init(size_t ts) {
    table_size = 1;
    while (table_size < ts) table_size <<= 1;
    if (table_size < 8) table_size = 8;
    data.assign(table_size * 3, 0);
    mask = table_size - 1;
  }
  
  /**
   * @brief Obtiene el índice (valor de 24 bits) almacenado en la posición lógica que se solicita
   * @param pos Posición lógica (bucket) en la tabla
   * @return El índice almacenado expandido a uint32_t
   */
  inline uint32_t get(size_t pos) const noexcept {
    size_t off = pos * 3;
    uint32_t v = (uint32_t)data[off] | ((uint32_t)data[off+1] << 8) | ((uint32_t)data[off+2] << 16);
    return v;
  }
  
  /**
   * @brief Almacena un índice de 24 bits en la posición lógica indicada
   * @param pos Posición lógica en la tabla
   * @param val Valor a empaquetar y almacenar. Debe de ser < 2^24
   */
  inline void set(size_t pos, uint32_t val) noexcept {
    size_t off = pos * 3;
    data[off]   = (uint8_t)(val & 0xFF);
    data[off+1] = (uint8_t)((val >> 8) & 0xFF);
    data[off+2] = (uint8_t)((val >> 16) & 0xFF);
  }
  
  /** @brief Vacía toda la tabla llenándola de ceros */
  inline void clear() noexcept { std::fill(data.begin(), data.end(), 0); }
};

// FlatSet ultra compacto, claves fijas: nwords uint64_t 

/**
 * @brief Implementación de un conjunto ultracompacto y plano basado en hashing lineal cerrado open addressing
 * @details Se evita los punteros de nodo y se reduce la fragmentación. Las claves son de tamaño fijo (nwords son palabras de 64 bits)
 */
struct FlatSetBitsCompact {
  size_t nwords;                 // Cantidad de palabras uint64_t que conforman cada clave
  vector<uint64> keys_packed;    // Almacenamiento contiguo de claves [key0_w0, key0_w1, key1_w0, key1_w1, ...]
  Packed24Table table;           // Tabla de resolución de colisiones que guarda los índices
  size_t count;                  // Cantidad de elementos actualmente almacenados
  
  FlatSetBitsCompact(): nwords(0), count(0) {}
  
  /**
   * @brief Inicializa el set con el tamaño de clave y una estimación de capacidad
   * @param _nwords Cantidad de uint64_t por clave
   * @param expected_elements Cantidad estimada de inserciones
   */
  void init(size_t _nwords, size_t expected_elements) {
    nwords = _nwords;
    count = 0;
    if (nwords == 0) throw runtime_error("nwords must be > 0");
    size_t need = expected_elements ? (size_t)ceil((double)expected_elements / MAX_LOAD_FACTOR) : 8;
    size_t table_size = 1;
    while (table_size < need) table_size <<= 1;
    if (table_size < 8) table_size = 8;
    table.init(table_size);
    // Reservar espacio para claves (nwords * expected_elements)
    if (expected_elements) keys_packed.reserve(expected_elements * nwords);
  }

  /** @brief Función mezcladora de 64 bits estándar */
  static inline uint64 mix_hash_uint64(uint64 x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return x;
  }

  /**
   * @brief Calcula el hash de una clave compuesta por varias palabras
   * @param words Puntero a la clave original
   * @return Valor hash final en 64 bits
   */
  uint64 hash_words(const uint64* words) const noexcept {
    uint64 h = 1469598103934665603ULL;
    for (size_t i = 0; i < nwords; ++i) {
      uint64 w = words[i];
      h ^= w + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    }
    return mix_hash_uint64(h);
  }

  /**
   * @brief Compara si dos claves representadas en memoria son iguales
   * @param words Puntero a la nueva clave a probar
   * @param idx Índice del elemento almacenado en keys_packed contra el que comparar
   * @return true si son idénticas palabra por palabra
   */
  bool key_equal_words(const uint64* words, size_t idx) const noexcept {
    const uint64* base = keys_packed.data() + idx * nwords;
    for (size_t i = 0; i < nwords; ++i) if (base[i] != words[i]) return false;
    return true;
  }

  /**
   * @brief Inserta un nuevo elemento si no existe, o retorna su índice si ya está
   * @param words Puntero a las palabras que definen la clave a insertar
   * @return pair donde el booleano es true si hubo inserción nueva, y size_t es el índice interno del registro
   */
  pair<bool,size_t> insert_or_find_from_words(const uint64* words) {
    // rehash si se supera la carga permitida
    if ((count + 1) > (size_t)(table.table_size * MAX_LOAD_FACTOR)) {
      rehash(table.table_size * 2);
    }
    uint64 h = hash_words(words);
    size_t pos = (size_t)h & table.mask;
    while (true) {
      uint32_t val = table.get(pos);
      if (val == 0) {
        // Añadir palabras al final
        size_t idx = count;
        keys_packed.insert(keys_packed.end(), words, words + nwords);
        ++count;
        // Almacenar idx+1 en tabla 24-bit (controlando overflow)
        if (idx + 1 >= (1u<<24)) throw runtime_error("index exceeds 24-bit capacity");
        table.set(pos, (uint32_t)(idx + 1));
        return {true, idx};
      } 
      else {
        size_t existing_idx = (size_t)(val - 1);
        if (existing_idx < count && key_equal_words(words, existing_idx)) {
          return {false, existing_idx};
        }
        pos = (pos + 1) & table.mask;
      }
    }
  }

  /**
   * @brief Aumenta el tamaño de la tabla hash y recoloca las entradas existentes
   * @param new_table_size Nuevo tamaño objetivo de la tabla. Esta debería de ser potencia de 2
   */
  void rehash(size_t new_table_size) {
    if (new_table_size < 8) new_table_size = 8;
    size_t ts = 1;
    while (ts < new_table_size) ts <<= 1;
    new_table_size = ts;
    Packed24Table newtab;
    newtab.init(new_table_size);
    size_t newmask = newtab.mask;
    for (size_t idx = 0; idx < count; ++idx) {
      const uint64* key = keys_packed.data() + idx * nwords;
      uint64 h = hash_words(key);
      size_t pos = (size_t)h & newmask;
      while (newtab.get(pos) != 0) pos = (pos + 1) & newmask;
      if (idx + 1 >= (1u<<24)) throw runtime_error("index exceeds 24-bit capacity during rehash");
      newtab.set(pos, (uint32_t)(idx + 1));
    }
    table.data.swap(newtab.data);
    table.table_size = newtab.table_size;
    table.mask = newtab.mask;
  }
};

/**
 * @brief Estructura de particionado para concurrencia sin bloqueos
 * @details Al aislar los estados generados en particiones diferentes con un hash, varios hilos pueden
 * trabajar al mismo tiempo añadiendo datos a sus respectivos sets sin bloquear constantemente un mutex global
 */
struct Shard {
  FlatSetBitsCompact fset;     // Conjunto hash exclusivo de este shard
  vector<uint32> frontier;     // Lista de índices en fset que conforman la frontera BFS actual
  vector<uint32> next_frontier;// Lista temporal para construir la frontera del siguiente paso
  std::mutex mtx;              // Mutex para proteger inserciones de hilos cruzados durante la fase de volcado
};

/**
 * @brief Función hash global para asignar un estado a su shard correspondiente
 * @param words Puntero al array de uint64_t del estado
 * @param nwords Cantidad de palabras
 * @return Un identificador hash global único
 */
uint64 hash_words_global(const uint64* words, size_t nwords) {
  uint64 h = 1469598103934665603ULL;
  for (size_t i = 0; i < nwords; ++i) {
    uint64 w = words[i];
    h ^= w + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  }
  h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
  return h;
}

/**
 * @brief Función principal que organiza el algoritmo BFS paralelo para descubrir todos los semigrupos dado un F
 * @details Se divide el espacio de búsqueda en shards asignados por hash. Se emplea buffers locales para que
 * la exploración se realice libre de bloqueos mutuos, aplicando los cambios bajo cerrojos por partición
 * @param F El número de Frobenius fijo
 * @throws runtime_error si el tamaño en palabras es excesivamente grande o F es inválido
 */
void generaSemigrupos_con_shards(int F) {
  if (F < 0) throw runtime_error("F must be non-negative");
  unsigned max_bit = 2U * (unsigned)F + 1U;
  size_t nwords = (size_t)(max_bit / 64U) + 1U;
  if (nwords == 0) nwords = 1;
  if (nwords > (1u<<20)) throw runtime_error("nwords excesivo (protección)");

  int num_threads = max(1, omp_get_max_threads());
  int num_shards = num_threads;
  vector<Shard> shards(num_shards);

  size_t per_shard_est = expected_states_main ? max((size_t)1, expected_states_main / (size_t)num_shards) : 1024;
  for (int s = 0; s < num_shards; ++s) shards[s].fset.init(nwords, per_shard_est);

  // Semigrupo inicial
  vector<uint16> S0 = semigrupoInicial(F);
  S0 = minimizarGeneradores(S0);
  vector<uint64> wtmp(nwords,0ULL);
  genvec_to_words(S0, wtmp.data(), nwords, max_bit);
  uint64 h0 = hash_words_global(wtmp.data(), nwords);
  int owner0 = (int)(h0 % (uint64)num_shards);
  auto pr0 = shards[owner0].fset.insert_or_find_from_words(wtmp.data());
  shards[owner0].frontier.push_back( (uint32)pr0.second );

  size_t count_totales = 1;
  vector<size_t> count_internos_sh(num_shards, 0);
  vector<size_t> count_hojas_sh(num_shards, 0);
  bool initTodoMenor = true;
  for (uint16 g : S0) if (g >= F) { initTodoMenor = false; break; }
  if (initTodoMenor) count_hojas_sh[owner0]++; else count_internos_sh[owner0]++;

  // Buffers por hilo. Vector de vectores planas de uint64 (cada clave ocupa nwords uint64)
  int max_threads = num_threads;
  // pre-allocate outer container
  vector<vector<uint64>> per_thread_owner_buf((size_t)max_threads * (size_t)num_shards);

  while (true) {
    bool any_frontier = false;
    for (int s = 0; s < num_shards; ++s) if (!shards[s].frontier.empty()) { any_frontier = true; break; }
    if (!any_frontier) break;

    // Limpiamos buffers
    for (auto &v : per_thread_owner_buf) { v.clear(); }

    #pragma omp parallel num_threads(max_threads)
    {
      int tid = omp_get_thread_num();
      int shard_id = tid % num_shards;
      Shard &sh = shards[shard_id];
      vector<uint64> tmpw(nwords);
      // Procesamos frontier local
      for (uint32 blk_idx : sh.frontier) {
        if ((size_t)blk_idx >= sh.fset.count) continue;
        const uint64* stored = sh.fset.keys_packed.data() + (size_t)blk_idx * nwords;
        // Encontrar primer bit
        int m_first = -1;
        for (size_t wi = 0; wi < nwords; ++wi) {
          uint64 w = stored[wi];
          if (w) {
            #if defined(__GNUG__) || defined(__clang__)
              int tz = __builtin_ctzll(w);
              m_first = (int)(wi*64 + tz);
            #else
              for (int b=0;b<64;++b) if ((w>>b)&1ULL) { m_first = (int)(wi*64 + b); break; }
            #endif
              break;
          }
        }
        int m_bound = (m_first >= 0) ? m_first : (F + 1);
        for (int x = 2; x < m_bound; ++x) {
          if (x == F) continue;
          size_t wi = (size_t)x >> 6;
          unsigned bi = x & 63;
          if (wi >= nwords) continue;
          if ((stored[wi] >> bi) & 1ULL) continue;
          // tmpw = stored | (1<<x)
          for (size_t i=0;i<nwords;++i) tmpw[i]=stored[i];
          tmpw[wi] |= (1ULL<<bi);
          // Convertimos a generadores para minimizar y validar el número de Frobenius
          vector<uint16> gens = words_to_genvec(tmpw.data(), nwords);
          gens = minimizarGeneradores(gens);
          if (!frobeniusValido(gens, F)) continue;
          // Recomputamos las palabras desde gens minimizado
          vector<uint64> finalw(nwords,0ULL);
          genvec_to_words(gens, finalw.data(), nwords, max_bit);
          uint64 h = hash_words_global(finalw.data(), nwords);
          int owner = (int)(h % (uint64)num_shards);
          // append finalw to per_thread_owner_buf[tid * num_shards + owner]
          auto &buf = per_thread_owner_buf[(size_t)tid * (size_t)num_shards + (size_t)owner];
          buf.insert(buf.end(), finalw.begin(), finalw.end());
        }
      }
    } // omp

    // Fase de inserción por propietario, con lock
    #pragma omp parallel for schedule(dynamic)
    for (int owner = 0; owner < num_shards; ++owner) {
      Shard &sh = shards[owner];
      std::lock_guard<std::mutex> lg(sh.mtx);
      for (int tid = 0; tid < max_threads; ++tid) {
        auto &flatbuf = per_thread_owner_buf[(size_t)tid * (size_t)num_shards + (size_t)owner];
        size_t off = 0, total = flatbuf.size();
        while (off + nwords <= total) {
          const uint64* key = flatbuf.data() + off;
          auto pr = sh.fset.insert_or_find_from_words(key);
          if (pr.first) {
            sh.next_frontier.push_back((uint32)pr.second);
            // Clasificación interno/hoja. Convertimos a gens solo para clasificación
            vector<uint16> gens = words_to_genvec(key, nwords);
            bool todoMenor = true;
            for (uint16 v : gens) if (v >= (uint16)F) { todoMenor = false; break; }
            if (todoMenor) count_hojas_sh[owner]++; else count_internos_sh[owner]++;
            #pragma omp atomic
            count_totales++;
          }
          off += nwords;
        }
        // Liberamos buffer
        vector<uint64>().swap(flatbuf);
      }
    }

    // swap frontiers
    for (int s = 0; s < num_shards; ++s) {
      shards[s].frontier.swap(shards[s].next_frontier);
      shards[s].next_frontier.clear();
      shards[s].next_frontier.shrink_to_fit();
    }
  } // while BFS

  size_t count_internos = 0, count_hojas = 0;
  for (int s = 0; s < num_shards; ++s) { count_internos += count_internos_sh[s]; count_hojas += count_hojas_sh[s]; }
  size_t computed_total = count_internos + count_hojas;
  cout << "\nTotal internos: " << count_internos
       << "   Total hojas: " << count_hojas << "\n";
  cout << "(Total semigrupos generados: " << computed_total << ")\n";
}

/**
 * @brief Función principal del programa
 * @details Se solicita al usuario el número de Frobenius F y se calcula el tiempo necesario para resolver la generación
 * @return Código de salida (0 si todo fue bien, 1 o 2 en errores)
 */
int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  cout << "Numero de Frobenius (F): ";
  string input;
  getline(cin, input);
  std::regex pattern("^[0-9]+$");
  if (!regex_match(input, pattern)) {
    cerr << "Entrada no valida.\n"; return 1;
  }
  int F = stoi(input);
  if (F <= 0) { cerr << "F debe ser positivo\n"; return 1; }

  // En caso de conocer una estimación de estados, la deberemos de poner para reservar memoria
  // expected_states_main = 800000000;
  omp_set_num_threads(omp_get_max_threads());
  auto t0 = chrono::high_resolution_clock::now();
  try {
    generaSemigrupos_con_shards(F);
  } 
  catch (const exception &e) {
    cerr << "Error: " << e.what() << "\n";
    return 2;
  }
  auto t1 = chrono::high_resolution_clock::now();
  auto dur = chrono::duration_cast<chrono::seconds>(t1 - t0).count();
  cout << "\nTiempo: " << dur << " s\n";
  return 0;
}

