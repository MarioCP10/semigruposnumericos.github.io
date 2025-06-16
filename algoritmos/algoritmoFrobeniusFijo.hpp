/**
 * @file algoritmoFrobeniusFijo.hpp
 * @brief Declaraciones para generar e imprimir semigrupos num�ricos con un valor fijo de Frobenius.
 * @author Mario Casas P�rez
 * @date 2025-06-16
 * @license MIT
 * @note Compilar con -std=c++11 o superior.
 */

#ifndef ALGORITMO_FROBENIUS_FIJO_HPP
#define ALGORITMO_FROBENIUS_FIJO_HPP

#include <string>
#include <vector>
#include <set>

namespace semigrupo {

/**
 * @brief Convierte un vector de generadores en su notaci�n de semigrupo <g1,�,gk>.
 * @param generadores Conjunto de generadores del semigrupo.
 * @return Cadena con el semigrupo en notaci�n angular.
 */
std::string semigrupoAString(const std::vector<int>& generadores);

/**
 * @brief Alias de semigrupoAString para contextos de Ap�ry.
 * @param generadores Conjunto de generadores.
 * @return Misma salida que semigrupoAString.
 */
std::string corchetesAngulares(const std::vector<int>& generadores);

/**
 * @brief Comprueba si un valor es representable como combinaci�n no negativa de generadores.
 * @param valor Entero a verificar.
 * @param generadores Conjunto de generadores del semigrupo.
 * @return true si existe combinaci�n de enteros >= 0 que sume `valor`.
 */
bool esRepresentable(int valor, const std::vector<int>& generadores);

/**
 * @brief Reduce un conjunto de generadores eliminando los redundantes.
 * @param generadores Vector inicial de generadores.
 * @return Nuevo vector sin generadores representables por los dem�s.
 */
std::vector<int> minimizarGeneradores(std::vector<int> generadores);

/**
 * @brief Verifica la validez del Frobenius de un semigrupo.
 * @param generadores Conjunto de generadores minimizados.
 * @param F Valor de Frobenius a comprobar.
 * @return true si F no es representable y F+1 s� lo es.
 */
bool frobeniusValido(const std::vector<int>& generadores, int F);

/**
 * @brief Construye el semigrupo inicial S0 = <F+1, F+2, �, 2F+1>.
 * @param F Valor de Frobenius.
 * @return Vector de generadores de S0.
 */
std::vector<int> semigrupoInicial(int F);

/**
 * @brief Candidato para a�adir nuevo generador.
 */
struct Candidato {
    int x;                         //Elemento candidato a a�adir.
    std::vector<int> semigrupo;   //Semigrupo resultante tras la minimizaci�n.
};

/**
 * @brief Obtiene los candidatos que pueden a�adirse manteniendo Frobenius = F.
 * @param S Conjunto actual de generadores minimizados.
 * @param F Valor de Frobenius objetivo.
 * @return Vector de Candidato con cada posible x y su semigrupo.
 */
std::vector<Candidato> obtenerCandidatos(const std::vector<int>& S, int F);

/**
 * @brief Calcula el conjunto de Ap�ry de un semigrupo respecto a un periodo.
 * @param generadores Conjunto de generadores minimizados.
 * @param periodo Entero de referencia (normalmente F+1).
 * @return Vector size=periodo con los menores valores representables de cada resto.
 */
std::vector<int> conjuntoApery(const std::vector<int>& generadores, int periodo);

/**
 * @brief Convierte el conjunto de Ap�ry a cadena.
 * @param generadores Conjunto de generadores.
 * @param periodo Periodo usado en el c�lculo (F+1).
 * @return Cadena con formato Ap(<generadores>, periodo) = {a0,�,a_{periodo-1}}.
 */
std::string aperyAString(const std::vector<int>& generadores, int periodo);

/**
 * @brief Comprueba si un semigrupo ya ha sido procesado para evitar duplicados.
 * @param S Conjunto de generadores minimizados.
 * @param vistos Conjunto de representaciones en cadena ya mostradas.
 * @return true si semigrupoAString(S) est� en vistos.
 */
bool yaVisto(const std::vector<int>& S, const std::set<std::string>& vistos);

/**
 * @brief Genera e imprime todos los semigrupos con Frobenius = F en BFS.
 * @param F Valor de Frobenius deseado.
 * @return Vector de todos los semigrupos minimizados encontrados.
 */
std::vector<std::vector<int>> generaSemigruposConF(int F);

}

#endif

