#include <stdio.h>

int main(){
	int ini;
	while(scanf("%d", &ini)!= EOF){
		
		printf("if (i > %d && i < %d)\n",ini, ini+600 );
		printf("  pulso(i) = pulso(i-1)+0.05;\n");
		printf("endif\n");
		
		printf("if (i >= %d && i < %d)\n", ini+600,ini+900);
		printf("  pulso(i) = pulso(i-1)-0.1;\n");
		printf("endif\n");
	}
	return 0;
}