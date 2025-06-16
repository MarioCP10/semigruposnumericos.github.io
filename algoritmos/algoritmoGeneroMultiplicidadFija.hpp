/**
 * @file algoritmoGeneroMultiplicidadFija.hpp
 * @brief Declaraciones para encontrar semigrupos num�ricos internos y hojas
 *        para un g�nero y multiplicidad dados.
 * @author Mario Casas P�rez
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
 * @brief Calcula el m�ximo com�n divisor entre dos enteros.
 * @param a Primer operando.
 * @param b Segundo operando.
 * @return El m�ximo com�n divisor de a y b.
 */
int maxCD(int a, int b);

/**
 * @brief Genera una clave ordenada para memoizaci�n a partir de un vector.
 * @param S Conjunto de generadores.
 * @return Cadena con los elementos de S ordenados y separados por comas.
 */
std::string clave(const std::vector<int>& S);

/**
 * @brief Calcula el conductor de un semigrupo num�rico.
 * @param S Conjunto de generadores.
 * @return El menor n tal que todos los enteros >= n son representables.
 */
int calculaConductor(const std::vector<int>& S);

/**
 * @brief Calcula el g�nero de un semigrupo num�rico.
 * @param S Conjunto de generadores.
 * @return N�mero de enteros no representables por S.
 */
int calculaGenero(const std::vector<int>& S);

/**
 * @brief Calcula el n�mero de Frobenius de un semigrupo generado por S.
 * @param S Conjunto de generadores.
 * @return M�ximo entero no representable (o -1 si todos son representables).
 */
int calculaFrobenius(const std::vector<int>& S);

/**
 * @brief Comprueba si S es minimal en forma expandida de Hilbert.
 * @param S Conjunto de generadores.
 * @return true si ning�n generador es redundante y se alcanza un bloque m�nimo.
 */
bool esMinimalHilbert(const std::vector<int>& S);

/**
 * @brief Determina si S genera un semigrupo hoja.
 * @param S Conjunto de generadores.
 * @param genero G�nero objetivo.
 * @return true si S genera el g�nero dado y su Frobenius supera al m�ximo generador.
 */
bool esHoja(const std::vector<int>& S, int genero);

/**
 * @brief Genera recursivamente todas las combinaciones de longitud fija.
 * @param subconjuntos Acumula las combinaciones generadas.
 * @param combinacion Vector auxiliar para la combinaci�n actual.
 * @param numeros Conjunto de n�meros candidatos.
 * @param comienzo �ndice desde el que a�adir en esta llamada.
 * @param tamano Longitud deseada de cada combinaci�n.
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
 * @param internos Lista de semigrupos num�ricos internos.
 * @param hojas    Lista de semigrupos num�ricos hoja.
 */
void comparaCantidades(
    const std::vector<std::vector<int>>& internos,
    const std::vector<std::vector<int>>& hojas
);

/**
 * @brief Busca y clasifica semigrupos num�ricos de g�nero y multiplicidad fijos.
 * @param genero        G�nero objetivo.
 * @param multiplicidad M�nimo elemento del semigrupo (multiplicidad).
 */
void encontrarSemigruposYHojas(int genero, int multiplicidad);

} 

#endif

