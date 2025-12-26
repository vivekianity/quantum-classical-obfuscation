from imports import *


class QCircuit:
    def __init__(self, *args, from_qasm=True):
        self.qc = args[0] if from_qasm else QuantumCircuit(*args)
        self.gate_mappings = {
            "X": self.qc.x,
            "Y": self.qc.y,
            "H": self.qc.h,
            "Z": self.qc.z,
            "S": self.qc.s,
            "T": self.qc.t,
            "I": self.qc.id,
            "CX": self.qc.cx,
            "CY": self.qc.cy,
            "CZ": self.qc.cz,
            "SWAP": self.qc.swap,
            "CSWAP": self.qc.cswap,
            "CCX": self.qc.ccx,
        }

        self.gate_index = {
            0: "X",
            1: "CX",
            2: "SWAP",
            3: "CCX",
            4: "CSWAP",
            5: "S",
        }

        self.qubit_count = {
            "X": 1,
            "CX": 2,
            "SWAP": 2,
            "CCX": 3,
            "CSWAP": 3,
            "S": 1,
        }

    @classmethod
    def from_qasm2(cls, file_name: str):
        qc = read_qasm2(file_name)
        return cls(qc)

    @classmethod
    def from_qasm3(cls, file_name: str):
        qc = read_qasm3(file_name)
        return cls(qc)

    def to_qasm2(self, file_name: str):
        write_qasm2(self.qc, file_name)

    def to_qasm3(self, file_name: str):
        write_qasm3(self.qc, file_name)

    def encrypt(self, key_size: int, effective_qubits: int = -1) -> str:
        rem = key_size
        key = []
        if effective_qubits == -1:
            effective_qubits = self.qc.num_qubits

        while rem > 0:
            idx = random.randint(0, len(list(self.gate_index.keys())) - 1)
            gate = self.gate_index[idx]
            contr = self.qubit_count[gate] + 1

            if (rem - contr == 0 or rem - contr > 1) and contr - 1 <= effective_qubits:
                encryptor: list[int | str] = ["@"]
                encryptor.append(idx)
                encryptor.append("#")

                qubits = random.sample(range(effective_qubits), contr - 1)
                self.gate_mappings[gate](*qubits)
                for qubit in qubits:
                    encryptor.append(qubit)
                    encryptor.append("|")

                if len(encryptor):
                    encryptor.pop()
                # encryptor.extend(qubits)
                key.append(encryptor)
                rem -= contr

        flattened_list = list(map(lambda x: str(x), chain.from_iterable(key)))
        return "".join(flattened_list[1:])

    def decrypt(self, key: str, incorrect_measure: str, effective_qubits: int = -1):
        deque = Deque()
        measurement = [int(x) for x in incorrect_measure]

        n = (
            effective_qubits if effective_qubits != -1 else len(measurement)
        )  # Use measurement length if unset

        if len(measurement) == 0:
            raise ValueError("Not a valid measurement result!")

        for item in key.split("@"):
            index, qubits = item.split("#")
            qubits = tuple(map(lambda x: int(x), qubits.split("|")))

            gate_idx = int(index)
            if gate_idx not in self.gate_index:
                raise ValueError(f"Invalid key {key}. Value at index {gate_idx+1} is invalid!")

            gate = self.gate_index[gate_idx]
            deque.push_back((gate, qubits))

        while not deque.empty():
            quantum_gate, qubits = deque.back()
            deque.pop_back()
            bit_pos = list(map(lambda x: n - x - 1, qubits))

            # Ensure bit_pos indices are within measurement bounds
            if all(0 <= pos < len(measurement) for pos in bit_pos):
                if quantum_gate == "CX":
                    ctrl, target = bit_pos
                    measurement[target] = measurement[target] ^ measurement[ctrl]

                elif quantum_gate == "CCX":
                    c1, c2, t = bit_pos
                    measurement[t] = (measurement[c1] & measurement[c2]) ^ measurement[
                        t
                    ]

                elif quantum_gate == "SWAP":
                    bit1, bit2 = bit_pos
                    measurement[bit2], measurement[bit1] = (
                        measurement[bit1],
                        measurement[bit2],
                    )

                elif quantum_gate == "CSWAP":
                    ctrl, bit1, bit2 = bit_pos
                    if measurement[ctrl] == 1:
                        measurement[bit2], measurement[bit1] = (
                            measurement[bit1],
                            measurement[bit2],
                        )

                elif quantum_gate == "X":
                    bit = bit_pos[0]
                    measurement[bit] ^= 1

        return "".join(list(map(str, measurement)))

    def measure(
        self,
        qubits_to_measure: list[int],
        add_bits: bool = False,
        measure_all: bool = True,
    ):
        if measure_all:
            self.qc.measure_all(add_bits=add_bits)
        else:
            if add_bits:
                # Create a new ClassicalRegister and add it to the circuit
                cr = ClassicalRegister(len(qubits_to_measure))
                self.qc.add_register(cr)
                # Map qubits to the new classical register
                self.qc.measure(
                    qubits_to_measure, [cr[i] for i in range(len(qubits_to_measure))]  # type: ignore
                )
            else:
                self.qc.measure(qubits_to_measure, qubits_to_measure)

    def draw(self, output="mpl", **kwargs):
        if output == "text":
            print(self.qc.draw(output=output, **kwargs))
        else:
            self.qc.draw(output=output, **kwargs)

    def compile_and_run(self, shots=1024):
        simulator = AerSimulator()
        transpiled = transpile(self.qc, simulator)
        result = simulator.run(transpiled, shots=shots).result()
        return transpiled, result.get_counts()

    @staticmethod
    def compute_tvd(original_res, obfus_res, shots):
        all_keys = set(original_res.keys()) | set(obfus_res.keys())
        return sum(
            abs(original_res.get(k, 0) - obfus_res.get(k, 0)) for k in all_keys
        ) / (2 * shots)

    @staticmethod
    def compute_dfc(obfus_res, correct_outputs, shots):
        preserved_count = sum(obfus_res.get(output, 0) for output in correct_outputs)
        max_wrong = max(
            (v for k, v in obfus_res.items() if k not in correct_outputs), default=0
        )
        return (preserved_count - max_wrong) / shots


def execute_grover(
    tvd_list,
    dfc_list,
    solution_sets: list[list[str]] | None = None,
    key_size: int = 20,
    iter: int = 1,
):
    def extract_marked_states(qasm_file_path: str) -> list[str]:
        with open(qasm_file_path, "r") as file:
            content = file.read()
        # The solutions are appended as a comment: "//Solutions = ['state1', 'state2', ...]"
        solutions_line = content.split("//Solutions = ")[-1].strip()
        # Parse the list of solutions (e.g., "['1010101010', '1100110011', ...]")

        marked_states = ast.literal_eval(solutions_line)
        return marked_states

    def compute_dfc(original: dict, encrypted: dict, marked_states: list[str]) -> float:
        total_shots = sum(encrypted.values())
        if total_shots == 0:
            return 1.0

        filtered_original = {
            state: count for state, count in original.items() if state in marked_states
        }

        correct_shots = 0
        for state, shots in encrypted.items():
            if state in filtered_original:
                correct_shots += shots

        incorrect_shots = {
            state: shots
            for state, shots in encrypted.items()
            if state not in filtered_original
        }

        incorrect_counts = sorted(incorrect_shots.values(), reverse=True)
        top_5_incorrect = (
            sum(incorrect_counts[:5])
            if len(incorrect_counts) >= 5
            else sum(incorrect_counts)
        )

        return (correct_shots - top_5_incorrect) / total_shots

    def compute_tvd(original: dict, encrypted: dict) -> float:
        total_shots = sum(original.values())
        if total_shots == 0:
            return 1.0

        all_states = set(original.keys()).union(set(encrypted.keys()))
        total_diff = 0
        for state in all_states:
            count_orig = original.get(state, 0)
            count_enc = encrypted.get(state, 0)
            total_diff += abs(count_orig - count_enc)

        return total_diff / (2 * total_shots)

    def generate_grover_circuits(solutions: list[list[str]], folder_path: str):
        grover = Grover()
        for idx in range(len(solutions)):
            qc = grover.create_grover_ciruit(solutions[idx])
            file_name = f"grover_{idx+1}"
            write_qasm3(qc, f"{folder_path}/{file_name}.qasm")
            with open(f"{folder_path}/{file_name}.qasm", "a") as file:
                file.write(f"\n\n//Solutions = {solutions[idx]}")

        print("Generated qasm files for the Grover circuits of provided solution sets.")

    qasm_folder_path = "qasm_files1/grover"

    if solution_sets is not None:
        if not os.path.exists(qasm_folder_path):
            os.makedirs(qasm_folder_path)
        generate_grover_circuits(solution_sets, qasm_folder_path)

    for file in os.listdir(qasm_folder_path):
        if os.path.isfile(os.path.join(qasm_folder_path, file)):
            qc = QCircuit.from_qasm3(os.path.join(qasm_folder_path, file))
            qc_copy = deepcopy(qc)

            qc_copy.measure(
                qubits_to_measure=list(range(qc_copy.qc.num_qubits)), add_bits=True
            )
            qc_copy.draw(
                filename=f"pics/grover/original_circuit_{iter}.png", output="mpl"
            )
            _, original_res = qc_copy.compile_and_run()
            plot_histogram(
                original_res,
                title="Original State (before encryption)",
                filename=f"pics/grover/original_result_{iter}.png",
            )
            print(f"Original result = {original_res}")

            key: str = qc.encrypt(key_size)

            qc.measure(
                qubits_to_measure=list(range(qc_copy.qc.num_qubits)), add_bits=True
            )
            qc.draw(
                filename=f"pics/grover/encrypted_circuit_{iter}.png",
                output="mpl",
            )
            _, encrypted_res = qc.compile_and_run()
            print(f"Incorrect result = {encrypted_res}")

            plot_histogram(
                encrypted_res,
                title="Incorrect matching state (after encryption)",
                filename=f"pics/grover/encrypted_result_{iter}.png",
            )

            corrected_res = {}
            for string, shots in encrypted_res.items():
                decrypted_measure = qc.decrypt(key, string)
                corrected_res[decrypted_measure] = shots

            plot_histogram(
                corrected_res,
                title="Actual measurement (decrypted)",
                filename=f"pics/grover/decrypted_result_{iter}.png",
            )
            marked_state = extract_marked_states(os.path.join(qasm_folder_path, file))
            dfc_list.append(compute_dfc(original_res, encrypted_res, marked_state))
            tvd_list.append(compute_tvd(original_res, encrypted_res))

    print(f"Length of tvd_list: {len(tvd_list)}")
    print(f"Length of dfc_list: {len(dfc_list)}")


def execute_bv(tvd_list, dfc_list, solution_sets: list[list[str]] | None = None):
    def generate_bv_circuits(solutions: list[list], folder_path: str):
        factor: str = ""
        bias: int = 0

        for item in solutions:
            factor = item[0]
            if len(item) == 2:
                bias = item[1]  # type: ignore

            bv = BV(factor, bias)
            for idx in range(len(solutions)):
                qc = bv.create_circuit(add_measurement=False)
                file_name = generate_unique_string()
                write_qasm2(qc, f"{folder_path}/{file_name}.qasm")
                with open(f"{folder_path}/{file_name}.qasm", "a") as file:
                    file.write(f"\n\n//Solutions = {solutions[idx]}")

        print(
            "Generated qasm files for the Bernstein Vazirani circuits of provided solution sets."
        )

    qasm_folder_path = "qasm_files/bv"

    if solution_sets is not None:
        if not os.path.exists(qasm_folder_path):
            os.makedirs(qasm_folder_path)
        generate_bv_circuits(solution_sets, qasm_folder_path)

    for file in os.listdir(qasm_folder_path):
        if os.path.isfile(os.path.join(qasm_folder_path, file)):
            qc = QCircuit.from_qasm2(os.path.join(qasm_folder_path, file))
            qc_copy = deepcopy(qc)

            qc_copy.measure(
                qubits_to_measure=list(range(qc.qc.num_qubits - 1)), measure_all=False
            )
            qc_copy.draw(filename="pics/bv/original_circuit.png")
            _, res_orig = qc_copy.compile_and_run(shots=1024)
            plot_histogram(
                res_orig,
                title="Actual marked state (without encryption)",
                color=["black", "white"],
            )
            plt.savefig("pics/bv/original_result.png")
            # plt.show()

            qc.qc.barrier()
            key: str = qc.encrypt(20, effective_qubits=qc.qc.num_qubits - 1)
            qc.qc.barrier()
            qc.measure(
                qubits_to_measure=list(range(qc.qc.num_qubits - 1)), measure_all=False
            )
            qc.draw(filename="pics/bv/encrypted_circuit.png", output="mpl", style="bw")
            _, encrypted_res = qc.compile_and_run(shots=20)
            plot_histogram(encrypted_res, title="Marked State after encryption")
            plt.savefig("pics/bv/encrypted_result.png")
            plt.show()
            corrected_res = {}

            for string, shots in encrypted_res.items():
                decrypted_measure = qc.decrypt(
                    key, string, effective_qubits=qc.qc.num_qubits - 1
                )
                corrected_res[decrypted_measure] = shots

            plot_histogram(corrected_res, title="Marked state after decryption")
            plt.savefig("pics/bv/decrypted_result.png")
            plt.show()
            tvd_list.append(tvd(original_res=res_orig, encrypted_res=encrypted_res))
            dfc_list.append(dfc(original_res=res_orig, encrypted_res=encrypted_res))


def execute_qaoa(
    tvd_list,
    dfc_list,
    adjacency_list: list[tuple[int, int, int | float]] | None = None,
    num_qubits: int = 5,
    reps: int = 5,
    idx: int = 1,
):
    plt.close("all")

    def generate_qaoa_circuit(
        adjacency_list: list[tuple[int, int, int | float]], folder_path: str
    ):
        qaoa = QAOA(num_qubits, adjacency_list)
        circuit = qaoa.create_circuit(reps=reps, measure_all=False)
        backend = AerSimulator()
        transpiled, result, objective_vals = qaoa.compile_and_run_circuit(
            circuit, backend, reps=reps, cost_hamiltonian=qaoa.cost_hamiltonian
        )
        file_name = f"{folder_path}/{generate_unique_string()}.qasm"
        write_qasm3(transpiled, file_name)
        with open(file_name, "a") as file:
            file.write(f"\n\n//Adjacency list = {adjacency_list}")

        print("Generated qasm files for the QAOA circuit of provided solution sets.")

    qasm_folder_path = "qasm_files/qaoa"

    if adjacency_list is not None:
        if not os.path.exists(qasm_folder_path):
            os.makedirs(qasm_folder_path)
        generate_qaoa_circuit(adjacency_list, qasm_folder_path)

    for file in os.listdir(qasm_folder_path):
        if os.path.isfile(os.path.join(qasm_folder_path, file)):
            qc = QCircuit.from_qasm3(os.path.join(qasm_folder_path, file))
            qc_copy = deepcopy(qc)

            SHOTS = 1024
            if idx == 1:
                qc.measure(
                    qubits_to_measure=list(range(qc.qc.num_qubits)), add_bits=True
                )
                qc.draw(
                    filename="pics/qaoa/original_circuit.png", output="mpl", style="bw"
                )
                global res
                _, res = qc.compile_and_run(shots=SHOTS)

                with open("original_result.txt", "a") as file:
                    file.write(f"{res}\n")

                plot_histogram(
                    res,
                    title="Original Result (without encryption)",
                    color=["black", "white"],
                )
                plt.savefig("pics/qaoa/original_result.png")

            key: str = qc_copy.encrypt(20)
            qc_copy.measure(
                qubits_to_measure=list(range(qc.qc.num_qubits)), add_bits=True
            )
            qc_copy.draw(
                filename=f"pics/qaoa/encrypted_circuit_{idx}.png",
                output="mpl",
                style="bw",
            )
            _, encrypted_res = qc_copy.compile_and_run(shots=SHOTS)

            with open("encrypted_res.txt", "a") as file:
                file.write(f"{idx}). {encrypted_res}\n")

            tvd = QCircuit.compute_tvd(res, encrypted_res, SHOTS)
            tvd_list.append(tvd)
            correct_res = [
                k
                for k, v in sorted(res.items(), key=lambda item: item[1], reverse=True)[
                    :4
                ]
            ]

            dfc = QCircuit.compute_dfc(encrypted_res, correct_res, SHOTS)

            dfc_list.append(dfc)

            plot_histogram(
                encrypted_res, title="Encrypted Result", color=["black", "white"]
            )
            # plt.show()
            plt.savefig(f"pics/qaoa/encrypted_result_{idx}.png")

            corrected_res = {}

            for string, shots in encrypted_res.items():
                decrypted_measure = qc.decrypt(key, string)
                corrected_res[decrypted_measure] = shots

            with open("decrypted_res.txt", "a") as file:
                file.write(f"{idx}). {corrected_res}\n")

            plot_histogram(
                corrected_res, title="Decrypted Result", color=["black", "white"]
            )
            # plt.show()
            plt.savefig(f"pics/qaoa/decrypted_result_{idx}.png")

    with open("tvd.txt", "w") as file:
        file.write(f"{tvd_list}")

    with open("dfc.txt", "w") as file:
        file.write(f"{dfc_list}")


def execute_hhl(
    hhl_tvd,
    hhl_dfc,
    initial_state: list[float] | None = None,
    clock_reg: int = 2,
    shots: int = 1024,
    generate_circuit: bool = False,
    iter: int = 1,
):
    plt.close("all")

    def generate_hhl_circuit(
        folder_path: str, initial_state: list[float] | None = None
    ):
        qc = HHL(clock_reg=clock_reg)

        if initial_state is not None:
            qc.qc.initialize(initial_state, qc.input_reg[0])
        else:
            qc.qc.x(qc.input_reg)

        qc.qc.h(qc.clock)
        qc.hhl()
        qc.qc.h(qc.clock)
        file_name = f"{folder_path}/hhl_{iter}.qasm"
        write_qasm3(qc.qc, file_name)
        with open(file_name, "a") as f:
            f.write(
                f"\n\n// Initial state = {initial_state if initial_state else '[0, 1]'}"
            )

        print("Generated QASM file for the HHL circuit.")

    qasm_folder_path = "qasm_files/hhl"

    if generate_circuit:
        if not os.path.exists(qasm_folder_path):
            os.makedirs(qasm_folder_path)
        generate_hhl_circuit(qasm_folder_path, initial_state)

    for file in os.listdir(qasm_folder_path):
        if os.path.isfile(os.path.join(qasm_folder_path, file)):
            qc = QCircuit.from_qasm3(os.path.join(qasm_folder_path, file))
            qc_copy = deepcopy(qc)

            qc.measure(
                qubits_to_measure=list(range(qc.qc.num_qubits)),
                measure_all=True,
                add_bits=True,
            )
            qc.draw(output="mpl", filename=f"pics/hhl/original_circuit_{iter}.png")
            _, res = qc.compile_and_run(shots=shots)

            with open(f"original_result_hhl_{iter}.txt", "a") as file:
                file.write(f"{res}\n")

            plot_histogram(res, title="Original Result (without encryption)")
            plt.savefig(f"pics/hhl/original_result_{iter}.png")
            plt.close("all")
            # plt.show()

            key = qc_copy.encrypt(key_size=20)
            qc_copy.measure(
                qubits_to_measure=list(range(qc_copy.qc.num_qubits)),
                measure_all=True,
                add_bits=True,
            )

            qc_copy.draw(
                output="mpl", filename=f"pics/hhl/encrypted_circuit_{iter}.png"
            )
            _, encrypted_res = qc_copy.compile_and_run(shots=shots)
            with open(f"encrypted_result_hhl_{iter}.txt", "a") as file:
                file.write(f"{encrypted_res}\n")

            plot_histogram(encrypted_res, title="Encrypted Result")
            plt.savefig(f"pics/hhl/encrypted_result_{iter}.png")
            plt.close()
            # plt.show()

            corrected_res = {}
            for string, shots in encrypted_res.items():
                decrypted_measure = qc.decrypt(key, string)
                corrected_res[decrypted_measure] = shots

            with open(f"decrypted_result_hhl_{iter}.txt", "a") as file:
                file.write(f"{corrected_res}\n")

            plot_histogram(corrected_res, title="Decrypted Result")
            plt.savefig(f"pics/hhl/decrypted_result_{iter}.png")
            plt.close()
            # plt.show()

            hhl_tvd.append(tvd(res, encrypted_res))
            hhl_dfc.append(dfc(res, encrypted_res))


def execute_shor(
    tvd_list,
    dfc_list,
    number_sets: list[tuple[int, int]] | None = None,
    key_size: int = 20,
    iter: int = 1,
):
    """Execute Shor's algorithm for given N and a pairs, saving QASM and results.

    Args:
        number_sets (list[tuple[int, int]]): List of (N, a) pairs to factorize
        key_size (int): Size of encryption key
    """

    def generate_shor_circuits(numbers: list[tuple[int, int]], folder_path: str):
        for idx, (N, a) in enumerate(numbers):
            shor_instance = Shor(N, a)
            qc = shor_instance.create_shor_circuit()
            file_name = f"shor_qasm_{iter}"
            write_qasm3(qc, f"{folder_path}/{file_name}.qasm")
            with open(f"{folder_path}/{file_name}.qasm", "a") as file:
                file.write(f"\n\n//N = {N}, a = {a}")
        print("Generated qasm files for the Shor circuits of provided number sets.")

    qasm_folder_path = "qasm_files/shor"

    if number_sets is not None:
        if not os.path.exists(qasm_folder_path):
            os.makedirs(qasm_folder_path)
        generate_shor_circuits(number_sets, qasm_folder_path)

    for file in os.listdir(qasm_folder_path):
        if os.path.isfile(os.path.join(qasm_folder_path, file)):
            file_path = os.path.join(qasm_folder_path, file)
            qc = QCircuit.from_qasm3(file_path)
            qc_copy = deepcopy(qc)

            with open(file_path, "r") as f:
                content = f.read()
            n_match = re.search(r"//N = (\d+), a = (\d+)", content)
            if n_match:
                N = int(n_match.group(1))
                a = int(n_match.group(2))
            else:
                raise ValueError(f"Could not extract N and a from {file_path}")
            shor_instance = Shor(N, a)

            # Original circuit execution
            qc_copy.qc.add_register(ClassicalRegister(shor_instance.n_count, "c"))
            qc_copy.measure(
                qubits_to_measure=list(range(shor_instance.n_count)), measure_all=False
            )
            qc_copy.draw(
                filename=f"pics/shor/original_circuit_{iter}.png",
                output="mpl",
                style="bw",
            )
            _, res_orig = qc_copy.compile_and_run()
            plot_histogram(
                res_orig,
                title="Original State (before encryption)",
                color=["black", "white"],
            )
            plt.savefig(f"pics/shor/original_result_{iter}.png")

            with open(f"results/shor/original_result_{iter}.txt", "a") as file:
                file.write(f"{iter}). {res_orig}\n")

            # Encrypted circuit execution
            key = qc.encrypt(key_size, effective_qubits=shor_instance.n_count)
            qc.qc.add_register(ClassicalRegister(shor_instance.n_count, "c"))
            qc.measure(
                qubits_to_measure=list(range(shor_instance.n_count)), measure_all=False
            )
            qc.draw(
                filename=f"pics/shor/encrypted_circuit_{iter}.png",
                output="mpl",
                style="bw",
            )
            _, encrypted_res = qc.compile_and_run()

            with open(f"results/shor/encrypted_result_{iter}.txt", "a") as file:
                file.write(f"{iter}). {encrypted_res}\n")

            tvd_val = tvd(res_orig, encrypted_res)
            dfc_val = dfc(res_orig, encrypted_res)

            if tvd_list is not None:
                tvd_list.append(tvd_val)

            if dfc_list is not None:
                dfc_list.append(dfc_val)

            plot_histogram(
                encrypted_res,
                title="Encrypted measurement (after encryption)",
                color=["black", "white"],
            )
            plt.savefig(f"pics/shor/encrypted_result_{iter}.png")
            plt.close("all")

            # Decrypted result
            corrected_res = {}
            for string, shots in encrypted_res.items():
                decrypted_measure = qc.decrypt(
                    key, string, effective_qubits=shor_instance.n_count
                )
                corrected_res[decrypted_measure] = shots
            plot_histogram(
                corrected_res, title="Decrypted measurement", color=["black", "white"]
            )
            plt.savefig(f"pics/shor/decrypted_result_{iter}.png")
            plt.close("all")

            with open(f"results/shor/decrypted_result_{iter}.txt", "a") as file:
                file.write(f"{iter}). {corrected_res}\n")

            # Post-processing
            measured_value = max(res_orig, key=res_orig.get)
            measured_int = int(measured_value, 2)
            phase = measured_int / (2**shor_instance.n_count)
            r = Fraction(phase).limit_denominator(shor_instance.N).denominator
            if r % 2 == 0:
                guesses = [
                    np.gcd(shor_instance.a ** (r // 2) - 1, shor_instance.N),
                    np.gcd(shor_instance.a ** (r // 2) + 1, shor_instance.N),
                ]
                for guess in guesses:
                    if guess != 1 and guess != shor_instance.N:
                        print(
                            f"Factors of {shor_instance.N} are: {guess} and {shor_instance.N // guess}"
                        )
                        break
            else:
                print(f"Found odd period r={r}; rerun may be needed")


if __name__ == "__main__":
    grover_solutions = [
        ["1011000110", "1110000011", "1011111111", "0111101001", "0010010111"],
        ["1101111110", "0111010010", "0000110101", "0100100010", "0100111100"],
        ["0111111111", "0111101000", "1100001101", "0011010111", "1000010000"],
        ["0010100000", "1101111100", "1001110101", "0010101100", "0110110001"],
        ["1011001101", "1010100001", "1011101100", "1010100111", "0100001110"],
        ["1110011100", "1000010001", "1101011010", "0010001110", "1010111100"],
        ["0011101000", "0010100110", "1011100010", "0100010111", "1001011101"],
        ["1100110100", "0100111101", "1100100001", "0100110011", "1011010010"],
        ["0101011000", "1000000010", "1011101010", "0010100101", "0110001100"],
        ["1001100001", "1010111011", "1010000100", "0111000010", "1010110101"],
    ]

    print("Executing Shor's algo")
    tvd_shor = []
    dfc_shor = []
    number_sets = [(15, 7)]

    execute_shor(tvd_shor, dfc_shor, number_sets)

    with open("tvd_shor.txt", "a") as file:
        file.write(f"{tvd_shor}\n")

    with open("dfc_shor.txt", "a") as file:
        file.write(f"{dfc_shor}\n")

    print("Shor's algo done")

    print("Executing Grover's algo")
    tvd_grover = []
    dfc_grover = []

    for i in range(1, 2):
        execute_grover(tvd_grover, dfc_grover, grover_solutions, iter=i)

    with open("tvd_grover.txt", "a") as file:
        file.write(f"{tvd_grover}")

    with open("dfc_grover.txt", "a") as file:
        file.write(f"{dfc_grover}")

    print("Grover's algo done")

    print("Executing BV algo")
    bv_tvd = []
    bv_dfc = []

    for i in range(1, 101):
        execute_bv(bv_tvd, bv_dfc)

    with open("bv_tvd.txt", "a") as file:
        file.write(f"{bv_tvd}\n")

    with open("bv_dfc.txt", "a") as file:
        file.write(f"{bv_dfc}\n")

    print("BV algo done")

    print("Executing QAOA algo")
    qaoa_tvd = []
    qaoa_dfc = []

    for i in range(1, 101):
        execute_qaoa(qaoa_tvd, qaoa_dfc, num_qubits=5, reps=2, idx=i)

    with open("tvd_qaoa.txt", "a") as file:
        file.write(f"{qaoa_tvd}\n")

    with open("dfc_qaoa.txt", "a") as file:
        file.write(f"{qaoa_dfc}\n")

    print("QAOA done")

    print("Executing HHL algo")
    hhl_tvd = []
    hhl_dfc = []
    execute_hhl(hhl_tvd, hhl_dfc, generate_circuit=True)

    with open("tvd_hhl.txt", "a") as file:
        file.write(f"{hhl_tvd}\n")

    with open("dfc_hhl.txt", "a") as file:
        file.write(f"{hhl_dfc}\n")

    print("HHL done")
