#include <eff.h>
#include <stdio.h>

int main()
{
    int i = 0;
    int set = 0;

    eff_pinmux_set(PINMUX_11, PINMUX_GPIO);
    eff_gpio_dir_set(GPIO_11, (GPIO_PIN_2 | GPIO_PIN_3), EFF_GPIO_OUT);
    eff_gpio_pull_set(GPIO_11, (GPIO_PIN_2 | GPIO_PIN_3), EFF_GPIO_PULL_NONE);

    while (1)
    {
        printf("Energy is everything! %i\r\n", i++);
        if (set)
        {
            eff_gpio_clear(GPIO_11, (GPIO_PIN_2 | GPIO_PIN_3));
            set = 0;
        }
        else
        {
            eff_gpio_set(GPIO_11, (GPIO_PIN_2 | GPIO_PIN_3));
            set = 1;
        }
        sleep(1);
    }
}
