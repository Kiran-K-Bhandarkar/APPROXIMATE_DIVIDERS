LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY REG IS
GENERIC
(
	INPUT_WIDTH: INTEGER := 16
);
PORT
(
	CLK, RST: IN STD_LOGIC;
	INPUT: IN STD_LOGIC_VECTOR(INPUT_WIDTH-1 DOWNTO 0);
	OUTPUT: OUT STD_LOGIC_VECTOR(INPUT_WIDTH-1 DOWNTO 0)
);
END ENTITY REG;

ARCHITECTURE BEHAVIORAL OF REG IS
BEGIN
  REG_PROCESS: PROCESS(CLK, RST, INPUT)
    BEGIN
	  IF RST = '1' THEN
	    OUTPUT <= (OTHERS => '0');
	  ELSIF RISING_EDGE(CLK) THEN
	    OUTPUT <= INPUT;
	  END IF;
  END PROCESS REG_PROCESS;
END ARCHITECTURE BEHAVIORAL;