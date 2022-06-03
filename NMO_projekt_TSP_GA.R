### NMO projekt:
#
### Problem - optymalizacja trasy kuriera rozwozacegocego paczki do wylosowanych paczkomatow Inpost na terenie Warszawy [TSP]
#
### Zastosowana metoda - algorytm genetyczny
#
### Grupa: 
#           Aleksandra Klos, 
#           Karolina Mierzwa, 
#           Jakub Malenczuk, 
#           Jakub Gorzala
#
### ___________________________________________________________________________________________________________________
  
  
  if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
  if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
  if (!require('arrow')) install.packages('arrow'); library('arrow')
  if (!require('geosphere')) install.packages('geosphere'); library('geosphere')
  if (!require('leaflet')) install.packages('leaflet'); library('leaflet')
  if (!require('sp')) install.packages('sp'); library('sp')
  if (!require('tibble')) install.packages('tibble'); library('tibble')
  if (!require('magrittr')) install.packages('magrittr'); library('magrittr')
  if (!require('GA')) install.packages('GA'); library('GA')
  if (!require('igraph')) install.packages('igraph'); library('igraph')
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  
  # Wczytywanie danych
  inpost <- read_parquet("punktyOdbioru.parquet", as_tibble = TRUE) %>%
              filter(brand == 'Inpost' & lat > 52.15 & lat < 52.35 & lng > 20.90 & lng < 21.15)
  inpost['id'] <- paste0('Inpost_box_', seq(nrow(inpost))) # Nadanie numeru id paczkomatom
  
  
  # Losowanie paczkomatow, ktore znajda sie na trasie - domyslna liczba = 20
  set.seed(142)
  n_boxes <- 20
  indices <- sample(nrow(inpost), n_boxes, replace = FALSE)
  post_box_delivery <- inpost[indices, ]
  
  
  # Mapa z wylosowanymi paczkomatami
  post_box_delivery_coordinates <- post_box_delivery
  coordinates(post_box_delivery_coordinates) <- ~lng+lat
  map_boxes <- leaflet(post_box_delivery_coordinates) %>% addMarkers(label = ~id, popup = ~id) %>% addTiles()
  map_boxes
  
  
  
  
  # Definiowanie funkcji _____________________________________________________________________________________
  
  
  initialisation <- function(inpost_df, n_pop){
    # Tworzenie populacji poczatkowej
    
    population_list <- vector(mode='list', length=n_pop)
    # Losowe poczatkowych rozwiazan
    for (i in 1:length(population_list)){
      population <- sample_n(inpost_df['id'], size=nrow(inpost_df), replace = FALSE)$id
      current_population <- list(order=population, distance_sum=NULL)
      population_list[[i]] <- current_population
    }
    return(population_list)
  }
  
  
  
  sum_distance_between_boxes <- function(inpost_df, population_list, distance_matrix){
    # Sumowanie odleglosci pomiedzy kolejnymi paczkomatami na trasie
    # dla wszystkich reprezentantow populacji
    
    for (population_index in 1:length(population_list)){
      current_population <- population_list[[population_index]]
      distance_current_population <- 0
      
      for (box_index in 1:length(current_population$order)){
        if (box_index != length(current_population$order)){
          first_box <- current_population$order[box_index]
          second_box <- current_population$order[box_index+1]
          distance_between_first_second_box <- distance_matrix[first_box,second_box]
          distance_current_population <- distance_current_population + distance_between_first_second_box
        }
      }
      population_list[[population_index]]['distance_sum'] <- distance_current_population
    }
    return(population_list)
  }
  
  
  
  selection_operator <- function(population_list, elite_size){
    # Selekcja metoda tarczy ruletki
    
    population_distances <- c()
    for (population_index in 1:length(population_list)){
      population_distances <- c(population_distances, population_list[[population_index]]$distance_sum)
    }
    population_distance_min <- min(population_distances)
    solution_fitness <- population_distance_min/population_distances
    solution_fitness <- solution_fitness / sum(solution_fitness)
    
    # Elita - wybor najlepszych reprezentantow do puli rodzicielskiej (pomijajac selekcje)
    population_fitness_sort <- sort(solution_fitness, decreasing = TRUE, index.return=TRUE)
    population_fitness_rank_elite <- population_fitness_sort$ix[1:elite_size]
    
    # Wybor pozostalych reprezentantow (przy wykorzystaniu operatora selekcji)
    left_population_indices <- setdiff(1:length(population_list), population_fitness_rank_elite)
    solution_fitness_cumsum <- cumsum(solution_fitness)
    selected_population_indices <- c()
    for(i in 1:length(left_population_indices)){
      selected_pop_index <- which.min(solution_fitness_cumsum[left_population_indices] < runif(1))
      selected_population_indices <- c(selected_population_indices, selected_pop_index)
    }  
    # Przypisanie wybranych reprezentantow do puli rodzicielskiej
    selected_population_list <- population_list[c(population_fitness_rank_elite, selected_population_indices)]
    
    return(selected_population_list)
  }
  
  
  
  crossover_operator_1 <- function(selected_population_list, elite_size){
    # Krzyzowanie - PMX
    
    new_population_list <- vector(mode='list', length=length(selected_population_list))
    
    # Pula rodzicielska (z wykluczeniem elity)
    pool_indices <- sample(setdiff(1:length(selected_population_list), 1:elite_size), 
                              size = (length(selected_population_list)-elite_size)/2, replace = FALSE)
    pool <- selected_population_list[pool_indices]
  
    if (elite_size != 0) {
    # Przypisywanie elity do nowej populacji (bez krzyzowania)
      for (elite_index in 1:elite_size){
        new_population_list[[elite_index]] <- selected_population_list[[elite_index]]['order']
      }
    }
    
    indices_odd <- seq(from=1, by=2, length.out=length(pool_indices))
    
    # Krzyzowanie pozostalych reprezentantow populacji
    for (population_index in 1:(length(pool_indices))){
      if (population_index != length(pool_indices)){
        parent_a <- pool[[population_index]]$order
        parent_b <- pool[[population_index+1]]$order
      }
      # Dwa punkty obciecia
      gene1_index <- ceiling(runif(n=1) * length(parent_a))
      gene2_index <- ceiling(runif(n=1) * length(parent_a))
      start_point <- min(max(min(gene1_index, gene2_index), 2), length(parent_a)-1)
      end_point <- max(min(max(gene1_index, gene2_index), length(parent_a)-1), 2)
  
      # Krok 1 - Swapping
      offspring_1 <- c(rep(0, start_point-1), parent_b[start_point:end_point], rep(0, length(parent_a)-(end_point)))
      offspring_2 <- c(rep(0, start_point-1), parent_a[start_point:end_point], rep(0, length(parent_a)-(end_point)))
  
      # Kork 2 - pozostawienie niezduplikowanych genow
      for (i in 1:length(offspring_1)){
        if ((parent_a[i] %!in% parent_b[start_point:end_point]) & (i %!in% seq(start_point, end_point))){
          offspring_1[i] <- parent_a[i]
        }
        if ((parent_b[i] %!in% parent_a[start_point:end_point]) & (i %!in% seq(start_point, end_point))){
          offspring_2[i] <- parent_b[i]
        }
      }
      
      # Krok 3 - mapowanie pozostalych genow
      mapping_matrix <- matrix(c(parent_a[start_point:end_point], parent_b[start_point:end_point]), nrow = 2, byrow = TRUE)
      for (ii in 1:length(offspring_1)){
        
        if (offspring_1[ii] == "0"){
          matrix_row_index <- 2
          matrix_column_index <- which(mapping_matrix[matrix_row_index,] == parent_a[ii])
          matrix_mapping_value <- mapping_matrix[matrix_row_index-1, matrix_column_index]
          if (matrix_mapping_value %!in% offspring_1){
            offspring_1[ii] <- matrix_mapping_value
          } else {
            while (matrix_mapping_value %in% offspring_1){
              matrix_column_index <- which(mapping_matrix[matrix_row_index,] == matrix_mapping_value)
              if (matrix_row_index == 2){matrix_row_index <- 1} else{matrix_row_index <- 2}
              matrix_mapping_value <- mapping_matrix[matrix_row_index, matrix_column_index]
              if (matrix_row_index == 2){matrix_row_index <- 1} else{matrix_row_index <- 2}
            }
            offspring_1[ii] <- matrix_mapping_value
          }
        }
        
        if (offspring_2[ii] == "0"){
          matrix_row_index <- 1
          matrix_column_index <- which(mapping_matrix[matrix_row_index,] == parent_b[ii])
          matrix_mapping_value <- mapping_matrix[matrix_row_index+1, matrix_column_index]
          if (matrix_mapping_value %!in% offspring_2){
            offspring_2[ii] <- matrix_mapping_value
          } else {
            while (matrix_mapping_value %in% offspring_2){
              matrix_column_index <- which(mapping_matrix[matrix_row_index,] == matrix_mapping_value)
              if (matrix_row_index == 2){matrix_row_index <- 1} else{matrix_row_index <- 2}
              matrix_mapping_value <- mapping_matrix[matrix_row_index, matrix_column_index]
              if (matrix_row_index == 2){matrix_row_index <- 1} else{matrix_row_index <- 2}
            }
            offspring_2[ii] <- matrix_mapping_value
          }
        }
      }
      
      # Przypisywanie offspringow do nowej populacji
      offspring_1_pop_index <- indices_odd[population_index]
      offspring_2_pop_index <- indices_odd[population_index] + 1
      new_population_list[[elite_size+offspring_1_pop_index]] <- list(order = offspring_1)
      new_population_list[[elite_size+offspring_2_pop_index]] <- list(order = offspring_2)
    }
    
    return(new_population_list)
  }
  
  
  
  crossover_operator_2 <- function(selected_population_list, elite_size){
    # Krzyzowanie v2
    
    new_population_list <- vector(mode='list', length=length(selected_population_list))
    # Pula rodzicielska (z wykluczeniem elity)
    pool_indices <- sample(setdiff(1:length(selected_population_list), 1:elite_size), size = length(selected_population_list)-elite_size, replace = FALSE)
    pool <- selected_population_list[pool_indices]
    # Przypisywanie elity do nowej populacji (bez krzyzowania)
    for (elite_index in 1:elite_size){
      new_population_list[[elite_index]] <- selected_population_list[[elite_index]]['order']
    }
    # Krzyzowanie pozostalych reprezentantow populacji
    for (population_index in 1:(length(pool_indices))){
      if (population_index != length(pool_indices)){
        parent_a <- pool[[population_index]]$order
        parent_b <- pool[[population_index+1]]$order
      } 
      gene1 <- ceiling(runif(n=1) * length(parent_a))
      gene2 <- ceiling(runif(n=1) * length(parent_a))
      start_point <- min(gene1, gene2)
      end_point <- max(gene1, gene2)
      
      # Tworzenie offspring'a
      sub_child_1 <- parent_a[start_point:end_point]
      sub_child_2 <- c()
      for (i in parent_b){
        if (i %!in% sub_child_1){
          sub_child_2 <- c(sub_child_2, i)
        }
      }
      offspring <- c(sub_child_1, sub_child_2)
      
      # Przypisywanie offspring'a do nowej populacji
      new_population_list[[population_index+elite_size]] <- list(order = offspring)
    }
    return(new_population_list)
  }
  
  
  
  mutation_operator <- function(new_population_list, mutation_rate){
    # Mutacja genow
    
    mutation_pop_n <- floor(length(new_population_list)*mutation_rate)
    for (population in new_population_list){
      for (q in 1:mutation_pop_n){
        gene_mut_a_index <- sample(1:length(population$order), 1)
        gene_mut_b_index <- sample(1:length(population$order), 1)
        
        gene_a <- population$order[gene_mut_a_index]
        gene_b <- population$order[gene_mut_b_index]
        
        population$order[gene_mut_a_index] <- gene_b
        population$order[gene_mut_b_index] <- gene_a
      }
    }
    
    return(new_population_list)
  }
  
  
  
  eval_population <- function(population_list){
    # Ewaluacja populacji - przypisuje najlepszy i sredni wynik 
    # dla reprezentantow danej populacji
    
    population_distances_vector <- c()
    for (population in population_list){
      distance_current_population <- population$distance_sum
      population_distances_vector <- c(population_distances_vector, distance_current_population)
    }
    min_distance_index <- which(population_distances_vector == min(population_distances_vector))[1]
    min_distance_population_order <- population_list[[min_distance_index]]$order
    
    return(list(min_distance=min(population_distances_vector),
                avg_distance=mean(population_distances_vector),
                min_distance_order=min_distance_population_order))
  }
  
  
  
  run_GA <- function(pop_size, elite_size, mutation_rate, n_iteration, post_box_delivery, population_list){
    # Uruchamianie algorytmu GA
    
    # Tworzenie macierzy zawierajacej odleglosci pomiedzy wszystkimi punktami
    post_box_coords <- post_box_delivery %>% select(lng, lat)
    distance_matrix <- as.matrix(
      distm(post_box_coords, fun = distHaversine)
    )
    rownames(distance_matrix) <- post_box_delivery$id
    colnames(distance_matrix) <- post_box_delivery$id
    
    # Losowanie populacji poczatkowej
    population_list <- initialisation(inpost_df=post_box_delivery, n_pop=pop_size)
    best_result_list <- vector(mode='list', length=n_iteration)
  
    # Uruchamianie iteracji algorytmu
    start_time <- Sys.time()
    for (iteration in 1:n_iteration){
      print(paste('Numer iteracji:', iteration))
      
      # Fitness
      population_list <- sum_distance_between_boxes(inpost_df=post_box_delivery, population_list=population_list, distance_matrix=distance_matrix)
      eval_list <- eval_population(population_list)
      best_result_list[[iteration]]$iter <- iteration
      best_result_list[[iteration]]$min_dist <- eval_list$min_distance
      best_result_list[[iteration]]$avg_dist <- eval_list$avg_distance
      best_result_list[[iteration]]$min_distance_order <- eval_list$min_distance_order
      
      cat(paste('min distance:', eval_list$min_distance), "\n")
      cat(paste('avg distance:', eval_list$avg_distance), "\n\n")
      
      # Operator selekcji
      selected_population_list <- selection_operator(population_list, elite_size=elite_size)
      # Operator krzyzowania
      crossover_list <- crossover_operator_1(selected_population_list, elite_size=elite_size)
      # Operator mutacji
      mutation_list <- mutation_operator(crossover_list, mutation_rate=mutation_rate)  
      # Przypisywanie nowej populacji
      population_list <- mutation_list
    }
    end_time <- Sys.time()
    print(paste('Execution total time:', end_time-start_time))
    
    return(best_result_list)
  }






# Uruchamianie algorytmu (wlasna implementacja) __________________________________________________________________________________________

# PARAMETRY
pop_size <- 150
elite_size <- 1
mutation_rate <- 0.2
n_iteration <- 500

# Uruchamianie algorytu GA
best_result_list <- run_GA(pop_size, elite_size, mutation_rate, n_iteration, post_box_delivery, population_list)

# Przeksztalcenie best_result_list na ramke danych 
results_df <- dplyr::bind_rows(best_result_list)
min_dist_df <- results_df[results_df['min_dist'] == min(results_df['min_dist']), ]
min_dist_df <- min_dist_df[min_dist_df['avg_dist'] == min(min_dist_df['avg_dist']), ]
min_dist <- min_dist_df[1,2]$min_dist
# Calkowity dystans dla najlepszego rozwiazania
min_dist

# Wykres przebiegu iteracji
plot_df <- results_df %>% select(iter, min_dist, avg_dist) %>% group_by_all() %>% distinct()
colors <- c("średni dystans w populacji" = "blue", "minimalny dystans w populacji" = "green")
GA_plot_1 <- ggplot(plot_df, aes(x=iter, y=min_dist, color='minimalny dystans w populacji')) + 
              geom_line() + geom_point(shape = 21) +
              geom_line(aes(x=iter, y=avg_dist, color='średni dystans w populacji'), linetype='dashed') + 
              geom_point(aes(x=iter, y=avg_dist,  color='średni dystans w populacji'), shape = 21) +
              ylab('Całkowity dystans w metrach\n') +
              labs(color='') +
              ggtitle('Przebieg iteracji algorytmu GA') +
              scale_color_manual(values=colors) +
              theme_minimal() +
              theme(plot.title = element_text(size=18, hjust=0.5))
GA_plot_1

# Mapa z rozwiazaniem
GA_solution_df <- min_dist_df
colnames(GA_solution_df) <- c('iter', 'min_dist', 'avg_dist', 'id')
GA_map_df <- GA_solution_df %>% left_join(post_box_delivery, by='id')
solution_map_1 <- leaflet(GA_map_df) %>% 
                  addMarkers(lng = ~lng, lat = ~lat, label = ~id, popup = ~id) %>% 
                  addTiles() %>% addPolylines(data = GA_map_df, lng = ~lng, lat = ~lat)
solution_map_1




# _______ GA TSP - rozwiazanie za pomoca pakietu GA _______________________________________

# Tworzenie macierzy zawierajacej odleglosci pomiedzy wszystkimi punktami
post_box_coords <- post_box_delivery %>% select(lng, lat)
distance_matrix <- as.matrix(distm(post_box_coords, fun = distHaversine))
rownames(distance_matrix) <- post_box_delivery$id
colnames(distance_matrix) <- post_box_delivery$id

# Distance Matrix
tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[, 2:1]
  sum(distMatrix[route])
}

# Obliczanie Fitness
tpsFitness <- function(tour, ...) 1/tourLength(tour, ...)

# Uruchamianie GA
GA.fit <- ga(type = "permutation", fitness = tpsFitness, distMatrix = distance_matrix, lower = 1, 
             upper = nrow(post_box_delivery), popSize = 10, maxiter = 1000, run = 200, pmutation = 0.05, 
             monitor = NULL)


GA_post_box_solution_indices <- as.vector(GA.fit@solution[1,])
GA_post_box_solution <- post_box_delivery[GA_post_box_solution_indices,]
GA_post_box_solution_distance <- tourLength(GA.fit@solution[1,], distance_matrix)
# Calkowity dystans dla najlepszego rozwiazania
GA_post_box_solution_distance

# Wykres przebiegu iteracji
plot(GA.fit)


# Mapa z rozwiazaniem
solution_map_2 <- leaflet(GA_post_box_solution) %>% 
                    addMarkers(lng = ~lng, lat = ~lat, label = ~id, popup = ~id) %>% 
                    addTiles() %>% 
                    addPolylines(data = GA_post_box_solution, lng = ~lng, lat = ~lat)
solution_map_2
