/**
 * @file algoritmoGeneroFijo.hpp
 * @brief Declaraciones para encontrar semigrupos num�ricos internos y hojas para un g�nero dado.
 * @author Mario Casas P�rez
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
 * @brief Calcula el m�ximo com�n divisor entre dos enteros.
 * @param a Primer operando.
 * @param b Segundo operando.
 * @return El m�ximo com�n divisor de a y b.
 */
int maxCD(int a, int b);

/**
 * @brief Genera una clave ordenada para memoizaci�n a partir de un vector de generadores.
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
 * @brief Comprueba si S define una hoja.
 * @param S Conjunto de generadores.
 * @param genero G�nero fijo.
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
 * @brief Genera todas las combinaciones de tama�o fijo a partir de un conjunto de n�meros.
 * @param subconjuntos Contenedor donde se almacenar�n las combinaciones resultantes.
 * @param combinacion Vector temporal que guarda la combinaci�n en construcci�n.
 * @param numeros Vector de entrada con los elementos de los que se extraen las combinaciones.
 * @param comienzo �ndice en `numeros` desde el cual se empieza a tomar elementos.
 * @param tamano Tama�o objetivo de cada combinaci�n.
 */
void generaCombinaciones(std::vector<std::vector<int>>& subconjuntos,
                          std::vector<int>& combinacion,
                          const std::vector<int>& numeros,
                          int comienzo,
                          int tamano);

/**
 * @brief Compara la cantidad de semigrupos num�ricos internos y hojas e imprime el resultado.
 * @param internos Vector de semigrupos num�ricos internos.
 * @param hojas   Vector de semigrupos num�ricos que son hojas.
 */
void comparaCantidades(const std::vector<std::vector<int>>& internos,
                       const std::vector<std::vector<int>>& hojas);

/**
 * @brief Encuentra e imprime todos los semigrupos internos y hojas para un g�nero fijo.
 * @param genero G�nero fijo dado.
 */
void encontrarSemigruposYHojas(int genero);

}

#endif

