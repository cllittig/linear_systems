#include "metodos/gauss_jordan.hpp"
#include <vector>
#include <stdexcept>
#include "matriz.hpp"



namespace gaussjordan
{
    std::vector<double> gauss_jordan(Matriz data, int current_row)
    {
        int row = data.getRows(), column = data.getColumns();


        //caso base
        if(row == 1 && column == 2){

            double response = data.getValue(current_row, row) / data.getValue(current_row, current_row) ;

            if(data.getValue(0, 0) == 0 && data.getValue(0,1) == 0)
                throw std::runtime_error("Sistema indeterminado");
            else if(data.getValue(0, 0) == 0 && data.getValue(0,1) != 0)
                throw std::runtime_error("Sistema impossivel");

            return {response};
        }


        //piv. parcial
        int max_row = current_row;
        for(int i= current_row+1; i< row;i++){
            if(abs(data.getValue(i, current_row )) > abs(data.getValue(max_row, current_row)))
                max_row=i;
        } 

        //troca. linhas
        if(max_row!= current_row){
            for(int j = current_row; j<=row; j++){
                data.trocarLinhas(current_row, max_row);
            }
        }

        //norm. linha pivo
        double pivo = data.getValue(current_row, current_row);
        for(int j = current_row; j<=row; j++){
            data.setValue(current_row, j, data.getValue(current_row,j)/pivo);
        }

        for(int i =0; i<row; i++){
            //não posso zerar a linha do pivo atual usando o pivo atual
            if(i==current_row) continue;

            double fator = data.getValue(i, current_row);

            //faz a normalização das linhas abaixo do atual pivo
            for(int k= current_row; k<=row; k++){
                double newValue = data.getValue(i, k) - fator * data.getValue(current_row, k);
                data.setValue(i,k, newValue);
            }
        }

        //chamada da recursão.
        std::vector<double> parcSolution = gauss_jordan(data, current_row+1);


        double currentSolution = data.getValue(current_row, row);

        std::vector<double> completeSolution;
        completeSolution.push_back(currentSolution);
        completeSolution.insert(completeSolution.end(),parcSolution.begin(), parcSolution.end());

        return completeSolution;
    }
};