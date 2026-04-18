/**
 * @file algoritmoGeneroFijo.hpp
 * @brief Declaraciones para encontrar semigrupos numéricos internos y hojas para un género dado.
 * @author Mario Casas Pérez
 * @date 2025-06-16
 * @license MIT
 * @note Compilar con -std=c++11 o superior.
 */

#ifndef ALGORITMO_GENERO_FIJO_HPP
#define ALGORITMO_GENERO_FIJO_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <set>

namespace semigrupo {

extern std::unordered_map<std::string, int> cacheGenero;
extern std::unordered_map<std::string, int> cacheFrobenius;

/**
 * @brief Calcula el máximo común divisor entre dos enteros.
 * @param a Primer operando.
 * @param b Segundo operando.
 * @return El máximo común divisor de a y b.
 */
int maxCD(int a, int b);

/**
 * @brief Genera una clave ordenada para memoización a partir de un vector de generadores.
 * @param S Conjunto de generadores.
 * @return Cadena con los elementos de S ordenados y separados por comas.
 */
std::string clave(const std::vector<int>& S);

/**
 * @brief Calcula el conductor de un semigrupo numérico.
 * @param S Conjunto de generadores.
 * @return El menor n tal que todos los enteros >= n son representables.
 */
int calculaConductor(const std::vector<int>& S);

/**
 * @brief Calcula el género de un semigrupo numérico.
 * @param S Conjunto de generadores.
 * @return Número de enteros no representables por S.
 */
int calculaGenero(const std::vector<int>& S);

/**
 * @brief Calcula el número de Frobenius de un semigrupo generado por S.
 * @param S Conjunto de generadores.
 * @return Máximo entero no representable (o -1 si todos son representables).
 */
int calculaFrobenius(const std::vector<int>& S);

/**
 * @brief Comprueba si S es minimal en forma expandida de Hilbert.
 * @param S Conjunto de generadores.
 * @return true si ningún generador es redundante y se alcanza un bloque mínimo.
 */
bool esMinimalHilbert(const std::vector<int>& S);

/**
 * @brief Comprueba si S define una hoja.
 * @param S Conjunto de generadores.
 * @param genero Género fijo.
 * @return true si S es hoja, false si no.
 */
bool esHoja(const std::vector<int>& S, int genero);

/**
 * @brief Comprueba si todos los pares en S son coprimos (mcd = 1).
 * @param S Conjunto de generadores.
 * @return true si mcd de cada par es 1.
 */
bool mcdEsUno(const std::vector<int>& S);

/**
 * @brief Genera todas las combinaciones de tamaño fijo a partir de un conjunto de números.
 * @param subconjuntos Contenedor donde se almacenarán las combinaciones resultantes.
 * @param combinacion Vector temporal que guarda la combinación en construcción.
 * @param numeros Vector de entrada con los elementos de los que se extraen las combinaciones.
 * @param comienzo Índice en `numeros` desde el cual se empieza a tomar elementos.
 * @param tamano Tamaño objetivo de cada combinación.
 */
void generaCombinaciones(std::vector<std::vector<int>>& subconjuntos,
                          std::vector<int>& combinacion,
                          const std::vector<int>& numeros,
                          int comienzo,
                          int tamano);

/**
 * @brief Compara la cantidad de semigrupos numéricos internos y hojas e imprime el resultado.
 * @param internos Vector de semigrupos numéricos internos.
 * @param hojas   Vector de semigrupos numéricos que son hojas.
 */
void comparaCantidades(const std::vector<std::vector<int>>& internos,
                       const std::vector<std::vector<int>>& hojas);

/**
 * @brief Encuentra e imprime todos los semigrupos internos y hojas para un género fijo.
 * @param genero Género fijo dado.
 */
void encontrarSemigruposYHojas(int genero);

}

#endif

