def create_subplot(data, ax=None):
    if ax is None:
        ax = plt.gca()
    more_data = do_something_on_data()  
    bp = ax.boxplot(more_data)
    return bp

# make figure with subplots
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))
create_subplot(data, ax1)

def plot_something(data, ax=None, **kwargs):
    ax = ax or plt.gca()
    # Do some cool data transformations...
    return ax.boxplot(data, **kwargs)
Then you can experiment with your plotting function by simply calling plot_something(my_data) and you can specify which axes to use like so.

fig, (ax1, ax2) = plt.subplots(2)
plot_something(data1, ax1, color='blue')
plot_something(data2, ax2, color='red')