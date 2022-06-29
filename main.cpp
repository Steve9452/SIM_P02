#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#define PAUSE int n; cin >> n;

using namespace std;

#include "data_structures/SDDS.h"
#include "geometry/mesh.h"
#include "gid/input_output.h"
#include "utilities/math_utilities.h"
#include "utilities/FEM_utilities.h"
template <typename A>
Data SDDS<A>::ref = Data();

void free_list(DS<DS<float>*>* L){
    //Se calcula la longitud de la lista
    int length;
    SDDS<DS<float>*>::extension(L, &length);
    //Se recorre la lista
    for(int i = 0; i < length; i++){
        //Se extrae la matriz actual
        DS<float>* temp;
        SDDS<DS<float>*>::extract(L,i,&temp);

        //Se libera el espacio en memoria de la matriz actual
        SDDS<float>::destroy(temp);
    }
    //Se libera el espacio en memoria de la lista de resultados
    SDDS<DS<float>*>::destroy(L);
}
// M[nnodes][nnodes] 
int main(int argc, char** argv){

    cout << "Initializing process...\nCreating auxiliar variables... ";
    // A = A inicial
    DS<float> *A, *T_full, *A_N, *M, *K1, *K2, *b;
    DS<int> *dirichlet_indices, *neumann_indices;

    DS<DS<float>*> *Result,*M_locals,*K1_locals,*b_locals, *K2_locals;

    //Se crea la lista para almacenar los resultados de cada tiempo
    SDDS<DS<float>*>::create(&Result,SINGLE_LINKED_LIST);

    cout << "OK\nReading input file and creating geometry object... ";

    //Se instancia un objeto Mesh
    Mesh* G = new Mesh();

    read_input_file(G, argv[1]);

    cout << "OK\nCreating temperature vectors... ";

    int nelems = G->get_quantity(NUM_ELEMENTS); //Se extrae la cantidad de elementos en la malla
    int nnodes = G->get_quantity(NUM_NODES);    //Se extrae la cantidad de nodos en la malla
  
    int free_nodes = nnodes - G->get_quantity(NUM_DIRICHLET_BCs);
    SDDS<float>::create(&A, free_nodes, 1, MATRIX);
    SDDS<float>::create(&T_full, nnodes, 1, MATRIX);
    SDDS<float>::create(&A_N, nnodes, 1, MATRIX);

    cout << "OK\nInitializing temperature vectors... ";
    // TODO: Cambiar initial temperature
    Math::init(A, G->get_parameter(INITIAL_WATER));

    SDDS<int>::create(&neumann_indices, G->get_quantity(NUM_NEUMANN_BCs), ARRAY);
    G->get_condition_indices(neumann_indices, NEUMANN);
    FEM::built_T_Neumann(A_N, G->get_parameter(NEUMANN_VALUE), neumann_indices);

    SDDS<int>::create(&dirichlet_indices, G->get_quantity(NUM_DIRICHLET_BCs), ARRAY);
    G->get_condition_indices(dirichlet_indices, DIRICHLET);
    float Ad = G->get_parameter(DIRICHLET_VALUE); //Se extrae el valor para las condiciones de Dirichlet
    FEM::build_full_T(T_full, A, Ad, dirichlet_indices);

    cout << "OK\nInitializing list of results... ";

    FEM::append_results(Result, T_full);

    cout << "OK\n\nObtaining time parameters and starting loop...\n";

    float dt = G->get_parameter(TIME_STEP);    //Se extrae el paso de tiempo
    float t = G->get_parameter(INITIAL_TIME);  //Se extrae el tiempo inicial
    t += dt;
    float tf = G->get_parameter(FINAL_TIME);   //Se extrae el tiempo final

    //Comienza el ciclo de ejecución, el cual continúa hasta alcanzar el tiempo final
    while( t <= tf ){

        cout << "\tWorking at TIME = " << t << "s:\n\n";
        
        //Se preparan los arreglos para almacenar todas las matrices locales de todos los elementos
        //La longitud de los 3 arreglos es igual a la cantidad de elementos
        SDDS<DS<float>*>::create(&M_locals, nelems, ARRAY);
        SDDS<DS<float>*>::create(&K1_locals, nelems, ARRAY);
        SDDS<DS<float>*>::create(&K2_locals, nelems, ARRAY);
        SDDS<DS<float>*>::create(&b_locals, nelems, ARRAY);

        //Se recorren los elementos
        for(int e = 0; e < nelems; e++){
            cout << "\t\tWorking with ELEMENT = " << e+1 << ":\n";
            //Se interpreta el contador como un ID de elemento, con la salvedad
            //que el contador comienza en 0 y los IDs comienzan en 1

            //Se extrae el elemento actual, se envía e+1 para compensar la diferencia en los conteos
            Element* current_elem = G->get_element(e+1);

            cout << "\t\tCalculating local systems... ";
            //Se calcula la M local y se añade al listado de matrices M. Se envían la densidad y el calor específico del material
            SDDS<DS<float>*>::insert(M_locals, e, FEM::calculate_local_M(current_elem));
            //Se calcula la K local y se añade al listado de matrices K. Se envía la conductividad térmica del material
            SDDS<DS<float>*>::insert(K1_locals, e, FEM::calculate_local_K_1(current_elem));
            //
            SDDS<DS<float>*>::insert(K2_locals, e, FEM::calculate_local_K_2(current_elem));
            //Se calcula la b local y se añade al listado de matrices b. Se envía la fuente de calor
            SDDS<DS<float>*>::insert(b_locals, e, FEM::calculate_local_b(current_elem));
            cout << "OK\n\n";
        }

        cout << "\tCreating global system...\n";

        //Se crean las matrices globales, y se inicializan todas sus posiciones con 0
        SDDS<float>::create(&M, nnodes, nnodes, MATRIX); Math::zeroes(M);
        SDDS<float>::create(&K1, nnodes, nnodes, MATRIX); Math::zeroes(K1);
        SDDS<float>::create(&K2, nnodes, nnodes, MATRIX); Math::zeroes(K2);
        SDDS<float>::create(&b, nnodes, 1, MATRIX);      Math::zeroes(b);

        //Se recorren los listados de matrices locales, un elemento a la vez
        for(int e = 0; e < nelems; e++){
            cout << "\t\tAssembling ELEMENT = " << e+1 << ":\n";
            //Se extrae el elemento actual, se envía e+1 para compensar la diferencia en los conteos
            Element* current_elem = G->get_element(e+1);
            DS<float> *temp;

            cout << "\t\tAssembling local matrices... ";
            //Se extrae la matriz M del elemento actual y se envía a ensamblaje
            SDDS<DS<float>*>::extract(M_locals,e,&temp);
            FEM::assembly(M, temp, current_elem, true);  //Se indica que ensamblará una matriz 3 x 3
            //Se extrae la matriz K1 del elemento actual y se envía a ensamblaje
            SDDS<DS<float>*>::extract(K1_locals,e,&temp);
            FEM::assembly(K1, temp, current_elem, true);
             //Se extrae la matriz K2 del elemento actual y se envía a ensamblaje
            SDDS<DS<float>*>::extract(K1_locals,e,&temp);
            FEM::assembly(K2, temp, current_elem, true);
            //Se extrae la matriz b del elemento actual y se envía a ensamblaje
            SDDS<DS<float>*>::extract(b_locals,e,&temp);
            FEM::assembly(b, temp, current_elem, false); //Se indica que ensamblará una matriz 3 x 1
            cout << "OK\n\n";
        }

        cout << "\tApplying Neumann conditions... ";
        //Se agrega la matriz de valores de Neumann a la matriz b global

        // b - A_N
        Math::sub_in_place(b,A_N);
        cout << "OK\n\n";

        cout << "\tApplying Dirichlet conditions.. ";
        //Se modifican las matrices globales para aplicar las condiciones de Dirichlet
        // FEM::apply_Dirichlet(nnodes, free_nodes, &b, K1, Ad, dirichlet_indices);
        // FEM::apply_Dirichlet(nnodes, free_nodes, &K1, dirichlet_indices);

        // FEM::apply_Dirichlet(nnodes, free_nodes, &b, K2, Ad, dirichlet_indices);
        // FEM::apply_Dirichlet(nnodes, free_nodes, &K2, dirichlet_indices);
        cout << "\t<<<<...>>>> ";
        FEM::apply_Dirichlet(nnodes, free_nodes, &M, dirichlet_indices);
        cout << "OK\n\n";

        cout << "\tCalculating temperature at next time step.\n\tUsing FEM generated formulas and Forward Euler... ";

        /*
            Se procede a ejecutar la ecuación de transferencia de calor en su versión discretizada con Forward Euler:

                        A^(i+1) = A^i + M^(-1) * delta_t * ( b - K * A^i )

                        A^(i+1) = A^i + M^(-1) * delta_t * ( b - K1 * A^i - K2 * A^i )

            En la expresión anterior, a la matriz b ya se le han incorporado el vector columna de las condiciones
            de Neumann, y el vector columna generado por la aplicación de las condiciones de Dirichlet.
        */

        //Se ejecuta K * A, donde A son las temperaturas en el tiempo actual
        DS<float>* temp = Math::product(K1,A);
        DS<float>* temp1 = Math::product(K2,A);
        //Se multiplica el contenido del resultado anterior por -1 para simular la resta
        Math::product_in_place(temp, -1);
        Math::product_in_place(temp1, -1);
        //Se suma el resultado de -K*A a la matriz b
        Math::sum_in_place(b, temp);
        Math::sum_in_place(b, temp1);
        //Se multiplica el contenido del resultado anterior por delta_t, el paso de tiempo
        Math::product_in_place(b, dt);
        //Se multiplica el resultado anterior por la inversa de la matriz M, y el resultado se añade a los
        //resultados del tiempo actual, obteniendo así los resultados del siguiente tiempo
        //SDDS<float>::show(M,false);
        //DS<float>* temp2 = Math::inverse(M);
        DS<float>* temp2 = Math::inverse_Cholesky(M);
        //SDDS<float>::show(temp2,false);
        DS<float>* temp3 = Math::product( temp2, b );
        Math::sum_in_place(A, temp3 );
        //La matriz temp ya no será utilizada, por lo que se libera su espacio en memoria
        SDDS<float>::destroy(temp);SDDS<float>::destroy(temp2);SDDS<float>::destroy(temp3);

        cout << "OK\n\n\tUpdating list of results... ";

        //Se construye la matriz de resultados completa para el tiempo actual
        FEM::build_full_T(T_full, A, Ad, dirichlet_indices);
        //Se añade la matriz de resultados complete del tiempo actual a la lista de resultados
        FEM::append_results(Result, T_full);

        cout << "OK\n\nCleaning up and advancing in time... ";

        //Se libera todo el espacio en memoria utilizado en el tiempo actual
        free_list(M_locals);
        free_list(K1_locals);
        free_list(K2_locals);
        free_list(b_locals);
        SDDS<float>::destroy(M);
        SDDS<float>::destroy(K1);
        SDDS<float>::destroy(K2);
        SDDS<float>::destroy(b);

        //Avanzamos al siguiente tiempo a calcular
        t = t + dt;

        cout << "OK\n\n";
    }

    cout << "Writing output file... ";

    //El ingreso en las listas enlazadas simples se hace al inicio en SDDS,
    //por lo que hay que invertir la lista para tener los resultados en el orden
    //en el que se calcularon
    SDDS<DS<float>*>::reverse(Result);
    //Se envía el listado de resultados, y el nombre de archivo recibido en línea de
    //comandos, para crear el archivo de salida que será utilizado por GiD
    write_output_file(Result, argv[1]);

    cout << "OK\n\nCleaning up and finalizing process... ";

    //Se libera el espacio en memoria asignado para todas las estructuras utilizadas
    SDDS<float>::destroy(A); SDDS<float>::destroy(T_full); SDDS<float>::destroy(A_N);
    SDDS<int>::destroy(dirichlet_indices); SDDS<int>::destroy(neumann_indices);
    //Para la lista de resultados, por estar compuesta por otras estructuras de datos,
    //es necesario liberar una por una:
    
    free_list(Result);

    //Se libera el objeto Mesh
    delete G;
    
    cout << "OK\n\nHave a nice day!! :D\n";
    
    return 0;
}
