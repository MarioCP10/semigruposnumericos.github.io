/**
 * @file algoritmoGeneroMultiplicidadFija.hpp
 * @brief Declaraciones para encontrar semigrupos numéricos internos y hojas
 *        para un género y multiplicidad dados.
 * @author Mario Casas Pérez
 * @date 2025-06-16
 * @license MIT
 * @note Compilar con -std=c++11 o superior.
 */

#ifndef ALGORITMO_GENERO_MULTIPLICIDAD_FIJA_HPP
#define ALGORITMO_GENERO_MULTIPLICIDAD_FIJA_HPP

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
 * @brief Genera una clave ordenada para memoización a partir de un vector.
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
 * @brief Determina si S genera un semigrupo hoja.
 * @param S Conjunto de generadores.
 * @param genero Género objetivo.
 * @return true si S genera el género dado y su Frobenius supera al máximo generador.
 */
bool esHoja(const std::vector<int>& S, int genero);

/**
 * @brief Genera recursivamente todas las combinaciones de longitud fija.
 * @param subconjuntos Acumula las combinaciones generadas.
 * @param combinacion Vector auxiliar para la combinación actual.
 * @param numeros Conjunto de números candidatos.
 * @param comienzo Índice desde el que añadir en esta llamada.
 * @param tamano Longitud deseada de cada combinación.
 */
void generaCombinaciones(
    std::vector<std::vector<int>>& subconjuntos,
    std::vector<int>& combinacion,
    const std::vector<int>& numeros,
    int comienzo,
    int tamano
);

/**
 * @brief Comprueba si los elementos de S son globalmente coprimos.
 * @param S Conjunto de generadores.
 * @return true si el mcd acumulado termina en 1.
 */
bool mcdEsUno(const std::vector<int>& S);

/**
 * @brief Compara la cantidad de semigrupos internos y hojas.
 * @param internos Lista de semigrupos numéricos internos.
 * @param hojas    Lista de semigrupos numéricos hoja.
 */
void comparaCantidades(
    const std::vector<std::vector<int>>& internos,
    const std::vector<std::vector<int>>& hojas
);

/**
 * @brief Busca y clasifica semigrupos numéricos de género y multiplicidad fijos.
 * @param genero        Género objetivo.
 * @param multiplicidad Mínimo elemento del semigrupo (multiplicidad).
 */
void encontrarSemigruposYHojas(int genero, int multiplicidad);

} 

#endif

