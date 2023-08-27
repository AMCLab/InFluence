                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 6 : Electron deposits sufficient energy and moves positively in x and y.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if ((x0 > x_dimension) && (y0 > y_dimension)) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));

                        // Check if the new coordinates are within the eh_charge_counter boundaries
                        if ((i_coordinate + translation_x <= perfect_image_0 - 1) && (j_coordinate + translation_y <= perfect_image_1 - 1)) {
                            eh_charge_counter.coeffRef(i_coordinate + translation_x, j_coordinate + translation_y) += new_eh_pairs;
                        }
                    }
                    
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 7 : Electron deposits sufficient energy and moves negatively in x and y.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if ((x0 < -x_dimension) && (y0 < -y_dimension)) {
                        // Electron moves negatively in x and y

                        

                        // Calculate the translations
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));

                        // Check if the new coordinates are within the eh_charge_counter boundaries
                        if ((i_coordinate - translation_x >= 0) && (j_coordinate - translation_y >= 0)) {
                            eh_charge_counter.coeffRef(i_coordinate - translation_x, j_coordinate - translation_y) += new_eh_pairs;
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 8 : Electron deposits sufficient energy and moves positively in x only.                        //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (x0 > x_dimension) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        if (i_coordinate + translation_x <= perfect_image_0 - 1) {
                            eh_charge_counter.coeffRef(i_coordinate + translation_x, j_coordinate) += new_eh_pairs;
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 9 : Electron deposits sufficient energy and moves negatively in x only.                        //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (x0 < -x_dimension) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        if (i_coordinate - translation_x >= 0) {
                            eh_charge_counter.coeffRef(i_coordinate - translation_x, j_coordinate) += new_eh_pairs;
                            
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 10 : Electron deposits sufficient energy and moves positively in y only.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (y0 > y_dimension) {
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));
                        if (j_coordinate + translation_y <= perfect_image_0 - 1) {
                            eh_charge_counter.coeffRef(i_coordinate, j_coordinate + translation_y) += new_eh_pairs;
                            
                           
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 11 : Electron deposits sufficient energy and moves negatively in y only.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (y0 < -y_dimension) {
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));
                        if (j_coordinate - translation_y >= 0) {
                            eh_charge_counter.coeffRef(i_coordinate, j_coordinate - translation_y) += new_eh_pairs;

                            
                        }
                   
                    }