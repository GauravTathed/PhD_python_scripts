    def evolve_state(self, verbose: bool = False):

        self.initialise_system()
 
        # define the value of the magnetic field change for this run

        self.per_run_mag_noise = self.generate_magnetic_noise()

        self.per_run_laser_noise = self.generate_laser_noise()

        self.per_run_freq_miscal = self.generate_frequency_miscalibration()

        self.per_run_pulse_miscal = self.generate_pulse_time_miscalibration()

        self.per_run_pulse_var = self.generate_pulse_time_variation()
 
        # for self.use_bussed_freq_miscalibrations, we take a single normal

        # distribution sample and multiply this by the uncertainty value for

        # each individual transition read from the .txt file

        self.bussed_freq_miscal_sample = np.random.normal()
 
        # we are ALWAYS using the gaussian sampled line noise contribution now

        # whether compensating the line signal or not

        self.per_run_line_noise = self.generate_gaussian_line_noise()
 
        # print("Noise values for this shot:")

        # print(f"  Mag: {self.per_run_mag_noise}")

        # print(f"  Laser: {self.per_run_laser_noise}")

        # print(f"  Freq miscal: {self.per_run_freq_miscal}")

        # print(f"  Pulse miscal: {self.per_run_pulse_miscal}")

        # print(f"  Pulse var: {self.per_run_pulse_var}")
 
        if self.line_triggering:

            self.exp_start_time = 0

        else:

            self.exp_start_time = np.random.uniform(0, 1e6 / 60)

        total_exp_time = self.exp_start_time
 
        if len(self.dict_evolution['transitions']) == 0:

            # without any transitions, just return identity

            total_unitary = np.eye(len(self.basis_order), dtype=complex)

        else:

            # build the overall unitary of the pulses

            unitaries = np.empty(

                (len(self.dict_evolution['transitions']), len(

                    self.basis_order), len(self.basis_order)),

                dtype=complex)
 
            # build the list of transition sensitivities in the matrix order to

            # be used for each unitary's diagonal for detuning error

            transitions_in_basis_order = np.array([[

                self.basis_order[0][1], self.basis_order[it][0],

                self.basis_order[it][1]

            ] for it in range(len(self.basis_order))])[1:]
 
            ordered_table_idx = self.convert_triplet_index(

                transitions_in_basis_order)

            ordered_trans_sens = np.array([

                self.transition_sensitivities[tuple(table_idx)]

                for table_idx in ordered_table_idx

            ])
 
            for idx, trans in enumerate(self.dict_evolution['transitions']):

                table_idx = self.convert_triplet_index(trans)

                angle = self.dict_evolution['angles'][

                    idx] / np.pi  # need angle first
 
                if self.dict_evolution['pulse_times'] is not None:

                    trans_pitime = self.dict_evolution['pulse_times'][idx] / (

                        angle)

                    inp_trans_pitime = trans_pitime

                else:

                    trans_pitime = self.transition_pitimes[tuple(

                        table_idx)]  # us

                    inp_trans_pitime = None
 
                # print('Transition', trans)

                # print('pi-time', trans_pitime)

                # breakpoint()
 
                trans_sens = self.transition_sensitivities[tuple(

                    table_idx)]  # MHz/G
 
                # set pulse time for this transition

                time = trans_pitime * angle
 
                # if verbose:

                #     if self.line_noise_compensated:

                #         print('Using Gaussian line noise compensation.')

                #     else:

                #         print(

                #             'Line noise integrand', 2 * np.pi * trans_sens *

                #             self.generate_line_noise_integrand(total_exp_time))

                #     print('trans sens', trans_sens)

                #     print('exp time', total_exp_time - self.exp_start_time)
 
                # first put together whether pulses are going from S->D or D->S

                phase = self.dict_evolution['phases'][idx]
 
                if self.dict_evolution['noise_flags'][idx] == 1:

                    # add in noise sources, starting with line signal (60/180

                    # Hz) if the line signal is compensated, we only have the

                    # contribution from the gaussian noise associated with the

                    # discrepancy between the fit line and the actual frequency

                    if self.line_noise_compensated:

                        phase += (2 * np.pi * trans_sens *

                                  self.per_run_line_noise *

                                  (total_exp_time - self.exp_start_time))
 
                        if idx == 0:

                            real_detuning = (

                                2 * np.pi * ordered_trans_sens *

                                self.generate_line_noise_detuning(

                                    (1e6 / 60) * np.random.random()))

                        else:

                            real_detuning = (2 * np.pi * ordered_trans_sens *

                                             self.per_run_line_noise)
 
                    else:

                        # when line signal is not compensated, we have the full

                        # contribution from the line signal - as well as the

                        # same gaussian noise but now added on top of the line

                        # signal contribution
 
                        # deterministic part

                        phase += (

                            2 * np.pi * trans_sens *

                            (self.generate_line_noise_integrand(total_exp_time)

                             - self.generate_line_noise_integrand(

                                 self.exp_start_time)))
 
                        # gaussian noisy part

                        phase += (2 * np.pi * trans_sens *

                                  self.per_run_line_noise *

                                  (total_exp_time - self.exp_start_time))
 
                        if idx == 0:

                            real_detuning = (

                                2 * np.pi * ordered_trans_sens *

                                self.generate_line_noise_detuning(

                                    (1e6 / 60) * np.random.random()))

                        else:

                            # deterministic part

                            real_detuning = (2 * np.pi * ordered_trans_sens *

                                             self.generate_line_noise_detuning(

                                                 total_exp_time))

                            # gaussian noisy part

                            real_detuning += (2 * np.pi * ordered_trans_sens *

                                              self.per_run_line_noise)
 
                    # print(2 * np.pi * trans_sens *

                    #     self.generate_line_noise_integrand(total_experiment_time))
 
                    # now the magnetic field Gaussian noise

                    real_detuning += (2 * np.pi * ordered_trans_sens *

                                      self.per_run_mag_noise)

                    phase += (2 * np.pi * trans_sens * self.per_run_mag_noise *

                              (total_exp_time - self.exp_start_time))
 
                    # then the laser noise (taken from the [0,2,0] decay fit)

                    real_detuning += np.full(

                        len(ordered_table_idx),

                        2 * np.pi * 1e-6 * self.per_run_laser_noise)

                    phase += (2 * np.pi * 1e-6 * self.per_run_laser_noise *

                              (total_exp_time - self.exp_start_time))
 
                    # next the frequency miscalibration this now involves a

                    # check that we either use bussed frequency estimations/not

                    if self.use_bussed_freq_calibrations:

                        trans_key = "[" + ",".join(map(str, trans)) + "]"

                        bussed_uncertainty_for_trans = self.bussed_freq_miscalibrations.get(

                            trans_key, {}).get('NewStd', None)

                        bussed_freq_miscal = (self.bussed_freq_miscal_sample *

                                              bussed_uncertainty_for_trans)

                        real_detuning += 2 * np.pi * 1e-6 * bussed_freq_miscal

                        phase += (2 * np.pi * 1e-6 * bussed_freq_miscal *

                                  (total_exp_time - self.exp_start_time))

                        # print('term bussed', 2 * np.pi * 1e-6 * bussed_freq_miscal_sample)

                        # print('term bussed phase', 2 * np.pi * 1e-6 *

                        #           bussed_freq_miscal_sample *

                        #           (total_exp_time - self.exp_start_time))
 
                        # real_detuning += 2 * np.pi * 1e-6 * self.per_run_freq_miscal

                        # phase += (2 * np.pi * 1e-6 * self.per_run_freq_miscal *

                        #           (total_exp_time - self.exp_start_time))
 
                        # print('term bussed', 2 * np.pi * 1e-6 * self.per_run_freq_miscal)

                        # print('term bussed phase', 2 * np.pi * 1e-6 *

                        #           self.per_run_freq_miscal *

                        #           (total_exp_time - self.exp_start_time))

                        # breakpoint()

                    else:

                        # next the frequency miscalibration

                        real_detuning += np.full(

                            len(ordered_table_idx),

                            2 * np.pi * 1e-6 * self.per_run_freq_miscal)

                        phase += (2 * np.pi * 1e-6 * self.per_run_freq_miscal *

                                  (total_exp_time - self.exp_start_time))
 
                    # print(

                    #     f"Real exp time: {total_exp_time-self.exp_start_time} us"

                    # )
 
                new_unit = self.build_unitary(

                    transition=trans,

                    time=time,

                    inp_trans_pitime=inp_trans_pitime,

                    phase=phase,

                    real_detuning=real_detuning,

                    pulse_miscal=self.per_run_pulse_miscal,

                    pulse_var=self.per_run_pulse_var)

                unitaries[idx, :, :] = new_unit
 
                # if idx == 0:

                #     pass

                # else:

                #     total_experiment_time += time

                total_exp_time += time
 
            total_unitary = unitaries[0, :, :]

            # print(self.transitions[0],self.angles[0], self.phases[0])

            for i in range(1, len(self.dict_evolution['transitions'])):

                # print(self.transitions[i],self.angles[i], self.phases[i])

                total_unitary = np.matmul(unitaries[i, :, :], total_unitary)
 
        # print("Simulator's unitary\n", total_unitary)

        output_vec = total_unitary @ self.init_vec

        # print(self.init_vec)

        # print(output_vec)
 
        final_state = np.square((np.abs(output_vec)))

        self.total_experiment_time = total_exp_time - self.exp_start_time

        # print('Final state', final_state)

        if verbose:

            print('Total experiment time:',

                  np.round(self.total_experiment_time, 2), 'us')
 
        return final_state

 