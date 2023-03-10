LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY EXACT_CELL_Q IS PORT
(
	X, Y, BIN: IN STD_LOGIC;
	D, BOUT	 : OUT STD_LOGIC
);
END ENTITY EXACT_CELL_Q;

ARCHITECTURE BEHAVIORAL OF EXACT_CELL_Q IS
BEGIN
	D 	 <=  X XOR Y XOR BIN;
	BOUT <= (BIN AND (NOT(X XOR Y))) OR ((NOT(X) AND Y));
END ARCHITECTURE BEHAVIORAL;

